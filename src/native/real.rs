use faddeeva_sys::*;

pub fn w_im(x: f64) -> f64 {
    unsafe { Faddeeva_w_im(x) }
}

pub fn erf(x: f64) -> f64 {
    unsafe { Faddeeva_erf1(x) }
}

pub fn erfi(x: f64) -> f64 {
    unsafe { Faddeeva_erfi1(x) }
}

pub fn erfc(x: f64) -> f64 {
    unsafe { Faddeeva_erfc1(x) }
}

pub fn erfcx(x: f64) -> f64 {
    unsafe { Faddeeva_erfcx1(x) }
}

pub fn dawson(x: f64) -> f64 {
    unsafe { Faddeeva_Dawson1(x) }
}
