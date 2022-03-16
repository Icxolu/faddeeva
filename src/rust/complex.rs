use super::ISPI;
use crate::rust::helper::{sinc, sinh_taylor, EXPA2N2};
use crate::Complex64;

pub fn w(Complex64 { re, im }: Complex64) -> Complex64 {
    if re == 0.0 {
        return Complex64 {
            re: super::real::erfcx(im),
            im: re,
        };
    } else if im == 0.0 {
        return Complex64 {
            re: (-(re * re)).exp(),
            im: super::real::w_im(re),
        };
    }
    let a = 0.518_321_480_430_086; // pi / sqrt(-log(eps*0.5))
    let c = 0.329_973_702_884_629_07; // (2/pi) * a;
    let a2 = 0.268_657_157_075_235_96; // a^2

    let x = re.abs();
    let y = im;
    let ya = im.abs();

    let ret; // return value
    let mut sum = [0.0; 5];

    let finish = |ret: Complex64, sum: [f64; 5]| {
        ret + Complex64 {
            re: 0.5 * c * y * (sum[1] + sum[2]),
            im: 0.5 * c * (sum[4] - sum[3]).copysign(re),
        }
    };

    // continued fraction is faster
    /* As pointed out by M. Zaghloul, the continued
    fraction seems to give a large relative error in
    Re w(z) for |x| ~ 6 and small |y|, so use
    algorithm 816 in this region: */
    if ya > 7.0 || (x > 6.0 && (ya > 0.1 || (x > 8.0 && ya > 1e-10) || x > 28.0)) {
        /* Poppe & Wijers suggest using a number of terms
            nu = 3 + 1442 / (26*rho + 77)
        where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
        (They only use this expansion for rho >= 1, but rho a little less
         than 1 seems okay too.)
        Instead, I did my own fit to a slightly different function
        that avoids the hypotenuse calculation, using NLopt to minimize
        the sum of the squares of the errors in nu with the constraint
        that the estimated nu be >= minimum nu to attain machine precision.
        I also separate the regions where nu == 2 and nu == 1. */
        let xs = if y < 0.0 { -re } else { re }; // compute for -z if y < 0
        if x + ya > 4000.0 {
            // nu <= 2
            if x + ya > 1e7 {
                // nu == 1, w(z) = i/sqrt(pi) / z
                // scale to avoid overflow
                if x > ya {
                    let yax = ya / xs;
                    let denom = ISPI / (xs + yax * ya);
                    ret = Complex64::new(denom * yax, denom);
                } else if ya.is_infinite() {
                    return if x.is_nan() || y < 0.0 {
                        Complex64::new(f64::NAN, f64::NAN)
                    } else {
                        Complex64::default()
                    };
                } else {
                    let xya = xs / ya;
                    let denom = ISPI / (xya * xs + ya);
                    ret = Complex64::new(denom, denom * xya);
                }
            } else {
                // nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5)
                let dr = xs * xs - ya * ya - 0.5;
                let di = 2.0 * xs * ya;
                let denom = ISPI / (dr * dr + di * di);
                ret = Complex64::new(denom * (xs * di - ya * dr), denom * (xs * dr + ya * di));
            }
        } else {
            // compute nu(z) estimate and do general continued fraction
            let (c0, c1, c2, c3, c4) = (3.9, 11.398, 0.08254, 0.1421, 0.2023); // fit
            let mut nu = (c0 + c1 / (c2 * x + c3 * ya + c4)).floor();
            let mut wr = xs;
            let mut wi = ya;

            nu = 0.5 * (nu - 1.0);
            loop {
                if nu <= 0.4 || nu.is_nan() {
                    break;
                }
                let denom = nu / (wr * wr + wi * wi);
                wr = xs - wr * denom;
                wi = ya + wi * denom;
                nu -= 0.5;
            }
            {
                // w(z) = i/sqrt(pi) / w:
                let denom = ISPI / (wr * wr + wi * wi);
                ret = Complex64::new(denom * wi, denom * wr);
            }
        }
        if y < 0.0 {
            // use w(z) = 2.0*exp(-z*z) - w(-z), but be careful of overflow in
            // exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
            // Note: rust handles complex exponential differently
            let c = if (ya - xs) * (xs + ya) == f64::NEG_INFINITY {
                Complex64::default()
            } else {
                Complex64::new((ya - xs) * (xs + ya), 2.0 * xs * y).exp()
            };
            2.0 * c - ret
        } else {
            ret
        }
    } else if x < 10.0 {
        let mut prod2ax = 1.0;
        let mut prodm2ax = 1.0;

        if y.is_nan() {
            return Complex64::new(y, y);
        }

        // use precomputed exp(-a2*(n*n)) table
        let expx2 = if x < 5e-4 {
            // compute sum4 and sum5 together as sum5-sum4
            let x2 = x * x;
            let expx2 = 1.0 - x2 * (1.0 - 0.5 * x2); // exp(-x*x) via Taylor
                                                     // compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
            let ax2 = 1.036_642_960_860_172 * x; // 2*a*x
            let exp2ax = 1.0 + ax2 * (1.0 + ax2 * (0.5 + 0.166_666_666_666_666_66 * ax2));
            let expm2ax = 1.0 - ax2 * (1.0 - ax2 * (0.5 - 0.166_666_666_666_666_66 * ax2));

            for n in 1.. {
                let coef = EXPA2N2[n - 1] * expx2 / (a2 * (n as f64 * n as f64) + y * y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum[0] += coef;
                sum[1] += coef * prodm2ax;
                sum[2] += coef * prod2ax;

                // really = sum5 - sum4
                sum[4] += coef * (2.0 * a) * n as f64 * sinh_taylor((2.0 * a) * n as f64 * x);

                // test convergence via sum3
                if coef * prod2ax < f64::EPSILON * sum[2] {
                    break;
                };
            }
            expx2
        } else {
            // x > 5e-4, compute sum4 and sum5 separately
            let expx2 = (-x * x).exp();
            let exp2ax = (2.0 * a * x).exp();
            let expm2ax = 1.0 / exp2ax;
            for n in 1.. {
                let coef = EXPA2N2[n - 1] * expx2 / (a2 * (n as f64 * n as f64) + y * y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum[0] += coef;
                sum[1] += coef * prodm2ax;
                sum[3] += (coef * prodm2ax) * (a * n as f64);
                sum[2] += coef * prod2ax;
                sum[4] += (coef * prod2ax) * (a * n as f64);
                // test convergence via sum5, since this sum has the slowest decay
                if (coef * prod2ax) * (a * n as f64) < f64::EPSILON * sum[4] {
                    break;
                };
            }
            expx2
        };

        // avoid spurious overflow for large negative y
        let expx2erfcxy = if y > -6.0 {
            // for y < -6, erfcx(y) = 2*exp(y*y) to double precision
            expx2 * super::real::erfcx(y)
        } else {
            2.0 * (y * y - x * x).exp()
        };
        if y > 5.0 {
            // imaginary terms cancel
            let sinxy = (x * y).sin();
            ret = Complex64::new(
                (expx2erfcxy - c * y * sum[0]) * (2.0 * x * y).cos()
                    + (c * x * expx2) * sinxy * sinc(x * y, sinxy),
                0.0,
            );
            finish(ret, sum)
        } else {
            let sinxy = (re * y).sin();
            let sin2xy = (2.0 * re * y).sin();
            let cos2xy = (2.0 * re * y).cos();
            let coef1 = expx2erfcxy - c * y * sum[0];
            let coef2 = c * re * expx2;
            ret = Complex64 {
                re: coef1 * cos2xy + coef2 * sinxy * sinc(re * y, sinxy),
                im: coef2 * sinc(2.0 * re * y, sin2xy) - coef1 * sin2xy,
            };
            finish(ret, sum)
        }
    } else {
        // x large: only sum3 & sum5 contribute (see above note)
        if x.is_nan() {
            return Complex64::new(x, x);
        }
        if y.is_nan() {
            return Complex64::new(y, y);
        }

        ret = Complex64 {
            re: (-x * x).exp(),
            im: 0.0,
        };

        // (round instead of ceil as in original paper; note that x/a > 1 here)
        let n0 = (x / a + 0.5).floor(); // sum in both directions, starting at n0
        let dx = a * n0 - x;
        sum[2] = (-dx * dx).exp() / (a2 * (n0 * n0) + y * y);
        sum[4] = a * n0 * sum[2];
        let exp1 = (4.0 * a * dx).exp();
        let mut exp1dn = 1.0;
        let mut dn = 1;
        while n0 - dn as f64 > 0.0 {
            // loop over n0-dn and n0+dn terms
            let np = n0 + dn as f64;
            let nm = n0 - dn as f64;
            let mut tp = (-(a * dn as f64 + dx) * (a * dn as f64 + dx)).exp();
            exp1dn *= exp1;
            let mut tm = tp * exp1dn; // trick to get tm from tp
            tp /= a2 * (np * np) + y * y;
            tm /= a2 * (nm * nm) + y * y;
            sum[2] += tp + tm;
            sum[4] += a * (np * tp + nm * tm);
            if a * (np * tp + nm * tm) < f64::EPSILON * sum[4] {
                return finish(ret, sum);
            }
            dn += 1;
        }
        loop {
            // loop over n0+dn terms only (since n0-dn <= 0)
            let np = n0 + dn as f64;
            dn += 1;
            let tp =
                (-(a * dn as f64 + dx) * (a * dn as f64 + dx)).exp() / (a2 * (np * np) + y * y);
            sum[2] += tp;
            sum[4] += a * np * tp;
            if a * np * tp < f64::EPSILON * sum[4] {
                return finish(ret, sum);
            }
        }
    }
}

pub fn erf(z @ Complex64 { re, im }: Complex64) -> Complex64 {
    if im == 0.0 {
        return Complex64 {
            re: super::real::erf(re),
            im, // preserve sign of 0
        };
    }

    if re == 0.0 {
        return Complex64 {
            re, // preserve sign of 0
            /* handle y -> Inf limit manually, since
            exp(y^2) -> Inf but Im[w(y)] -> 0, so
            IEEE will give us a NaN when it should be Inf */
            im: if im * im > 720.0 {
                if im > 0.0 {
                    f64::INFINITY
                } else {
                    f64::NEG_INFINITY
                }
            } else {
                (im * im).exp() * super::real::w_im(im)
            },
        };
    }

    let m_re_z2 = (im - re) * (re + im); // Re(-z^2), being careful of overflow
    let m_im_z2 = -2.0 * re * im; // Im(-z^2)

    if m_re_z2 < -750.0 {
        // underflow
        return if re >= 0.0 {
            Complex64 { re: 1.0, im: 0.0 }
        } else {
            Complex64 { re: -1.0, im: 0.0 }
        };
    }

    // Use Taylor series for small |z|, to avoid cancellation inaccuracy
    // erf(z) = 2/sqrt(pi) * z * (1 - z^2/3 + z^4/10 - z^6/42 + z^8/216 + ...)
    let taylor = || {
        let mz2 = Complex64::new(m_re_z2, m_im_z2); // -z^2
        z * (std::f64::consts::FRAC_2_SQRT_PI
            + mz2
                * (0.376_126_389_031_837_54
                    + mz2
                        * (0.112_837_916_709_551_26
                            + mz2 * (0.026_866_170_645_131_252 + mz2 * 0.005_223_977_625_442_188))))
    };

    /* for small |x| and small |xy|,
       use Taylor series to avoid cancellation inaccuracy:
         erf(x+iy) = erf(iy)
            + 2*exp(y^2)/sqrt(pi) *
              [ x * (1 - x^2 * (1+2y^2)/3 + x^4 * (3+12y^2+4y^4)/30 + ...
                - i * x^2 * y * (1 - x^2 * (3+2y^2)/6 + ...) ]
       where:
          erf(iy) = exp(y^2) * Im[w(y)]
    */
    let taylor_erfi = || {
        let (x2, y2) = (re * re, im * im);
        let expy2 = y2.exp();
        Complex64 {
            re: expy2
                * re
                * (std::f64::consts::FRAC_2_SQRT_PI
                    - x2 * (0.376_126_389_031_837_54 + 0.752_252_778_063_675_1 * y2)
                    + x2 * x2
                        * (0.112_837_916_709_551_26
                            + y2 * (0.451_351_666_838_205_05 + 0.150_450_555_612_735_02 * y2))),
            im: expy2
                * (super::real::w_im(im)
                    - x2 * im
                        * (std::f64::consts::FRAC_2_SQRT_PI
                            - x2 * (0.564_189_583_547_756_3 + 0.376_126_389_031_837_54 * y2))),
        }
    };

    /* Handle positive and negative x via different formulas,
    using the mirror symmetries of w, to avoid overflow/underflow
    problems from multiplying exponentially large and small quantities. */
    if re >= 0.0 {
        if re < 8e-2 {
            if im.abs() < 1e-2 {
                return taylor();
            } else if m_im_z2.abs() < 5e-3 && re < 5e-3 {
                return taylor_erfi();
            }
        }
        /* don't use complex exp function, since that will produce spurious NaN
        values when multiplying w in an overflow situation. */
        1.0 - m_re_z2.exp()
            * Complex64::from_polar(1.0, m_im_z2) /* {
                re: m_im_z2.cos(),
                im: m_im_z2.sin(),
            } */
            * w(Complex64::new(-im, re))
    } else {
        // x < 0
        if re > -8e-2 {
            // duplicate from above to avoid fabs(x) call
            if im.abs() < 1e-2 {
                return taylor();
            } else if m_im_z2.abs() < 5e-3 && re > -5e-3 {
                return taylor_erfi();
            }
        } else if re.is_nan() {
            return Complex64 {
                re: f64::NAN,
                im: if im == 0.0 { 0.0 } else { f64::NAN },
            };
        }
        /* don't use complex exp function, since that will produce spurious NaN
        values when multiplying w in an overflow situation. */

        m_re_z2.exp()
            * (Complex64 {
                re: m_im_z2.cos(),
                im: m_im_z2.sin(),
            } * w(Complex64::new(im, -re)))
            - 1.0
    }
}

/// erfi(z) = -i erf(iz)
pub fn erfi(Complex64 { re, im }: Complex64) -> Complex64 {
    let Complex64 { re, im } = erf(Complex64 { re: -im, im: re });
    Complex64 { re: im, im: -re }
}

/// erfc(z) = 1 - erf(z)
pub fn erfc(Complex64 { re, im }: Complex64) -> Complex64 {
    if re == 0.0 {
        return Complex64 {
            re: 1.0,
            /* handle y -> Inf limit manually, since
            exp(y^2) -> Inf but Im[w(y)] -> 0, so
            IEEE will give us a NaN when it should be Inf */
            im: if im * im > 720.0 {
                if im > 0.0 {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                }
            } else {
                -(im * im).exp() * super::real::w_im(im)
            },
        };
    }

    if im == 0.0 {
        return if re * re > 750.0 {
            // underflow
            Complex64 {
                re: if re >= 0.0 { 0.0 } else { 2.0 },
                im: -im, // preserve sign of 0
            }
        } else {
            Complex64 {
                re: if re >= 0.0 {
                    (-re * re).exp() * super::real::erfcx(re)
                } else {
                    2.0 - (-re * re).exp() * super::real::erfcx(-re)
                },
                im: -im, // preserve sign of 0
            }
        };
    }

    let m_re_z2 = (im - re) * (re + im);
    let m_im_z2 = -2.0 * re * im;

    if m_re_z2 < -750.0 {
        // underflow
        return if re >= 0.0 {
            Complex64::new(0.0, 0.0)
        } else {
            Complex64::new(2.0, 0.0)
        };
    }

    if re >= 0.0 {
        Complex64::new(m_re_z2, m_im_z2).exp() * w(Complex64::new(-im, re))
    } else {
        2.0 - Complex64::new(m_re_z2, m_im_z2).exp() * w(Complex64::new(im, -re))
    }
}

pub fn erfcx(Complex64 { re, im }: Complex64) -> Complex64 {
    w(Complex64 { re: -im, im: re })
}

pub fn dawson(z @ Complex64 { re: x, im: y }: Complex64) -> Complex64 {
    // handle axes separately for speed & proper handling of x or y = Inf or NaN
    if y == 0.0 {
        return Complex64 {
            re: super::real::w_im(x) / std::f64::consts::FRAC_2_SQRT_PI,
            im: -y, // preserve sign of 0
        };
    }

    if x == 0.0 {
        let y2 = y * y;
        return if y2 < 2.5e-5 {
            Complex64 {
                re: x, // preserve sign of 0
                im: y * (1.0 + y2 * (0.666_666_666_666_666_6 + y2 * 0.266_666_666_666_666_66)),
            }
        } else {
            Complex64 {
                re: x, // preserve sign of 0
                im: if y >= 0.0 {
                    y2.exp() - super::real::erfcx(y)
                } else {
                    super::real::erfcx(-y) - y2.exp()
                } / std::f64::consts::FRAC_2_SQRT_PI,
            }
        };
    }

    let m_re_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    let m_im_z2 = -2.0 * x * y; // Im(-z^2)
    let mz2 = Complex64::new(m_re_z2, m_im_z2); // -z^2

    // Use Taylor series for small |z|, to avoid cancellation inaccuracy
    // dawson(z) = z - 2/3 z^3 + 4/15 z^5 + ...
    let taylor = || z * (1.0 + mz2 * (0.666_666_666_666_666_6 + mz2 * 0.266_666_666_666_666_66));

    /*
    for small |y| and small |xy|,
    use Taylor series to avoid cancellation inaccuracy:
      dawson(x + iy)
       = D + y^2 (D + x - 2Dx^2)
           + y^4 (D/2 + 5x/6 - 2Dx^2 - x^3/3 + 2Dx^4/3)
       + iy [ (1-2Dx) + 2/3 y^2 (1 - 3Dx - x^2 + 2Dx^3)
             + y^4/15 (4 - 15Dx - 9x^2 + 20Dx^3 + 2x^4 - 4Dx^5) ] + ...
    where D = dawson(x)

    However, for large |x|, 2Dx -> 1 which gives cancellation problems in
    this series (many of the leading terms cancel).  So, for large |x|,
    we need to substitute a continued-fraction expansion for D.

       dawson(x) = 0.5 / (x-0.5/(x-1/(x-1.5/(x-2/(x-2.5/(x...))))))

    The 6 terms shown here seems to be the minimum needed to be
    accurate as soon as the simpler Taylor expansion above starts
    breaking down.  Using this 6-term expansion, factoring out the
    denominator, and simplifying with Maple, we obtain:

     Re dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / x
       = 33 - 28x^2 + 4x^4 + y^2 (18 - 4x^2) + 4 y^4
     Im dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / y
       = -15 + 24x^2 - 4x^4 + 2/3 y^2 (6x^2 - 15) - 4 y^4

    Finally, for |x| > 5e7, we can use a simpler 1-term continued-fraction
    expansion for the real part, and a 2-term expansion for the imaginary
    part.  (This avoids overflow problems for huge |x|.)  This yields:

    Re dawson(x + iy) = [1 + y^2 (1 + y^2/2 - (xy)^2/3)] / (2x)
    Im dawson(x + iy) = y [ -1 - 2/3 y^2 + y^4/15 (2x^2 - 4) ] / (2x^2 - 1)
    */
    let taylor_realaxis = || {
        let x2 = x * x;
        if x2 > 1600.0 {
            // |x| > 40
            let y2 = y * y;
            if x2 > 25e14 {
                // |x| > 5e7
                let xy2 = (x * y) * (x * y);
                return Complex64 {
                    re: (0.5 + y2 * (0.5 + 0.25 * y2 - 0.166_666_666_666_666_66 * xy2)) / x,
                    im: y
                        * (-1.0
                            + y2 * (-0.666_666_666_666_666_6 + 0.133_333_333_333_333_33 * xy2
                                - 0.266_666_666_666_666_66 * y2))
                        / (2.0 * x2 - 1.0),
                };
            }
            (1.0 / (-15.0 + x2 * (90.0 + x2 * (-60.0 + 8.0 * x2))))
                * Complex64 {
                    re: x * (33.0 + x2 * (-28.0 + 4.0 * x2) + y2 * (18.0 - 4.0 * x2 + 4.0 * y2)),
                    im: y * (-15.0 + x2 * (24.0 - 4.0 * x2) + y2 * (4.0 * x2 - 10.0 - 4.0 * y2)),
                }
        } else {
            let d = super::real::w_im(x) / std::f64::consts::FRAC_2_SQRT_PI;
            let y2 = y * y;
            Complex64 {
                re: d
                    + y2 * (d + x - 2.0 * d * x2)
                    + y2 * y2
                        * (d * (0.5 - x2 * (2.0 - 0.666_666_666_666_666_6 * x2))
                            + x * (0.833_333_333_333_333_4 - 0.333_333_333_333_333_3 * x2)),
                im: y
                    * (1.0 - 2.0 * d * x
                        + y2 * 0.666_666_666_666_666_6 * (1.0 - x2 - d * x * (3.0 - 2.0 * x2))
                        + y2 * y2
                            * (0.266_666_666_666_666_66
                                - x2 * (0.6 - 0.133_333_333_333_333_33 * x2)
                                - d * x
                                    * (1.0
                                        - x2 * (1.333_333_333_333_333_3
                                            - 0.266_666_666_666_666_66 * x2)))),
            }
        }
    };

    /* Handle positive and negative x via different formulas,
    using the mirror symmetries of w, to avoid overflow/underflow
    problems from multiplying exponentially large and small quantities. */
    if y >= 0.0 {
        if y < 5e-3 {
            if x.abs() < 5e-3 {
                return taylor();
            } else if m_im_z2.abs() < 5e-3 {
                return taylor_realaxis();
            }
        }
        let res = mz2.exp() - w(z);
        Complex64::new(-res.im, res.re) / std::f64::consts::FRAC_2_SQRT_PI
    } else {
        // y < 0
        if y > -5e-3 {
            // duplicate from above to avoid fabs(x) call
            if x.abs() < 5e-3 {
                return taylor();
            } else if m_im_z2.abs() < 5e-3 {
                return taylor_realaxis();
            }
        } else if y.is_nan() {
            return Complex64 {
                re: if x == 0.0 { 0.0 } else { f64::NAN },
                im: f64::NAN,
            };
        }
        let res = w(-z) - mz2.exp();
        Complex64::new(-res.im, res.re) / std::f64::consts::FRAC_2_SQRT_PI
    }
}
