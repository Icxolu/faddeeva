use crate::Complex64;
use faddeeva_sys::*;

pub fn w(z: Complex64) -> Complex64 {
    unsafe { Faddeeva_w(z.into(), 0.0) }.into()
}

pub fn erf(z: Complex64) -> Complex64 {
    unsafe { Faddeeva_erf(z.into(), 0.0) }.into()
}

pub fn erfi(z: Complex64) -> Complex64 {
    unsafe { Faddeeva_erfi(z.into(), 0.0) }.into()
}

pub fn erfc(z: Complex64) -> Complex64 {
    unsafe { Faddeeva_erfc(z.into(), 0.0) }.into()
}

pub fn erfcx(z: Complex64) -> Complex64 {
    unsafe { Faddeeva_erfcx(z.into(), 0.0) }.into()
}

pub fn dawson(z: Complex64) -> Complex64 {
    unsafe { Faddeeva_Dawson(z.into(), 0.0) }.into()
}
