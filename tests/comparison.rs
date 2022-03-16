#![cfg(feature = "native")]
mod common;
use common::{rel_err, C};

#[test]
fn w() {
    for &(z, _) in common::W {
        let w = faddeeva::complex::w(z);
        let fw = faddeeva::native::complex::w(z);
        let re_err = rel_err(w.re, fw.re);
        let im_err = rel_err(w.im, fw.im);
        assert!(re_err <= 1e-13);
        assert!(im_err <= 1e-13);
    }
}

#[test]
fn erf() {
    for &(z, _) in common::ERF {
        let w = faddeeva::complex::erf(z);
        let fw = faddeeva::native::complex::erf(z);
        let re_err = rel_err(w.re, fw.re);
        let im_err = rel_err(w.im, fw.im);

        assert!(re_err <= 1e-13);
        assert!(im_err <= 1e-13);
    }
}

#[test]
fn erfi() {
    // since erfi just calls through to erf, just one test should
    // be sufficient to make sure I didn't screw up the signs or something
    let z = C!(1.234, 0.5678);
    let w = faddeeva::complex::erfi(z);
    let fw = faddeeva::native::complex::erfi(z);

    let re_err = rel_err(w.re, fw.re);
    let im_err = rel_err(w.im, fw.im);

    assert!(re_err <= 1e-13);
    assert!(im_err <= 1e-13);
}

#[test]
fn erfcx() {
    // since erfcx just calls through to w, just one test should
    // be sufficient to make sure I didn't screw up the signs or something
    let z = C!(1.234, 0.5678);
    let w = faddeeva::complex::erfcx(z);
    let fw = faddeeva::native::complex::erfcx(z);

    let re_err = rel_err(w.re, fw.re);
    let im_err = rel_err(w.im, fw.im);

    assert!(re_err <= 1e-13);
    assert!(im_err <= 1e-13);
}

#[test]
fn erfc() {
    for &(z, _) in common::ERFC {
        let w = faddeeva::complex::erfc(z);
        let fw = faddeeva::native::complex::erfc(z);
        let re_err = rel_err(w.re, fw.re);
        let im_err = rel_err(w.im, fw.im);

        assert!(re_err <= 1e-13);
        assert!(im_err <= 1e-13);
    }
}

#[test]
fn dawson() {
    for &(z, _) in common::DAWSON {
        let w = faddeeva::complex::dawson(z);
        let fw = faddeeva::native::complex::dawson(z);
        let re_err = rel_err(w.re, fw.re);
        let im_err = rel_err(w.im, fw.im);

        assert!(re_err <= 1e-13);
        assert!(im_err <= 1e-13);
    }
}
