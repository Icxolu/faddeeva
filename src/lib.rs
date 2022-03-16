//! Computes various error functions ([erf](crate::complex::erf),
//! [erfc](crate::complex::erfc), [erfi](crate::complex::erfi),
//! [erfcx](crate::complex::erfcx)), including the
//! [Dawson](crate::complex::dawson) integral, in the complex plane, based on
//! algorithms for the computation of the [Faddeeva](crate::complex::w) function
//! `w(z) = exp(-z^2) * erfc(-i*z)`. Given `w(z)`, the error functions are
//! mostly straightforward to compute, except for certain regions where we have
//! to switch to Taylor expansions to avoid cancellation errors (e.g. near the
//! origin for `erf(z)`).
//!
//! To compute the Faddeeva function, we use a combination of two algorithms:
//!
//! For sufficiently large `|z|`, we use a continued-fraction expansion
//! for `w(z)` similar to those described in:
//!
//! > Walter Gautschi, "Efficient computation of the complex error function,"
//! > SIAM J. Numer. Anal. 7(1), pp. 187-198 (1970)
//!
//! > G. P. M. Poppe and C. M. J. Wijers, "More efficient computation of the
//! > complex error function," ACM Trans. Math. Soft. 16(1), pp. 38-46 (1990).
//!
//! Unlike those papers, however, we switch to a completely different algorithm
//! for smaller `|z|`:
//!
//! > Mofreh R. Zaghloul and Ahmed N. Ali, "Algorithm 916: Computing the
//! > Faddeyeva and Voigt Functions," ACM Trans. Math. Soft. 38(2), 15 (2011).
//!
//! # Note
//! This is a translation of Steven G. Johnson's
//! [Faddeeva](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) C/C++
//! library.
//!
//! # Features
//! * `native`: binds to the original C implementation and exposes the functions
//! under the `native` module. Mainly for testing purposes.

/* Copyright (c) 2012 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#[cfg(feature = "native")]
pub mod native;
mod rust;

pub use num_complex::Complex64;
pub use rust::*;
