#![cfg(feature = "native")]
mod common;
use common::{rel_err, C};

#[test]
fn w() {
    for &(z, w) in common::W {
        let fw = faddeeva::native::complex::w(z);
        let re_err = rel_err(w.re, fw.re);
        let im_err = rel_err(w.im, fw.im);

        assert!(re_err <= 1e-13);
        assert!(im_err <= 1e-13);
    }
}

#[test]
fn erf() {
    for &(z, w) in common::ERF {
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
    // erfi(z), computed with Maple
    let w = C!(1.081_032_284_405_373_2, 1.926_775_520_840_916_7);

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
    // erfcx(z), computed with Maple
    let w = C!(0.338_218_747_979_997_24, -0.111_607_747_081_164_85);

    let fw = faddeeva::native::complex::erfcx(z);
    let re_err = rel_err(w.re, fw.re);
    let im_err = rel_err(w.im, fw.im);

    assert!(re_err <= 1e-13);
    assert!(im_err <= 1e-13);
}

#[test]
fn erfc() {
    for &(z, w) in common::ERFC {
        let fw = faddeeva::native::complex::erfc(z);
        let re_err = rel_err(w.re, fw.re);
        let im_err = rel_err(w.im, fw.im);

        assert!(re_err <= 1e-13);
        assert!(im_err <= 1e-13);
    }
}

#[test]
fn dawson() {
    for &(z, w) in common::DAWSON {
        let fw = faddeeva::native::complex::dawson(z);
        let re_err = rel_err(w.re, fw.re);
        let im_err = rel_err(w.im, fw.im);

        assert!(re_err <= 1e-13);
        assert!(im_err <= 1e-13);
    }
}
