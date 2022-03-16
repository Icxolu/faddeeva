# Faddeeva

<!-- cargo-rdme start -->

Computes various error functions (erf,
erfc, erfi,
erfcx), including the
Dawson integral, in the complex plane, based on
algorithms for the computation of the Faddeeva function
`w(z) = exp(-z^2) * erfc(-i*z)`. Given `w(z)`, the error functions are
mostly straightforward to compute, except for certain regions where we have
to switch to Taylor expansions to avoid cancellation errors (e.g. near the
origin for `erf(z)`).

To compute the Faddeeva function, we use a combination of two algorithms:

For sufficiently large `|z|`, we use a continued-fraction expansion
for `w(z)` similar to those described in:

> Walter Gautschi, "Efficient computation of the complex error function,"
> SIAM J. Numer. Anal. 7(1), pp. 187-198 (1970)

> G. P. M. Poppe and C. M. J. Wijers, "More efficient computation of the
> complex error function," ACM Trans. Math. Soft. 16(1), pp. 38-46 (1990).

Unlike those papers, however, we switch to a completely different algorithm
for smaller `|z|`:

> Mofreh R. Zaghloul and Ahmed N. Ali, "Algorithm 916: Computing the
> Faddeyeva and Voigt Functions," ACM Trans. Math. Soft. 38(2), 15 (2011).

# Note
This is a translation of Steven G. Johnson's
[Faddeeva](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) C/C++
library.

# Features
* `native`: binds to the original C implementation and exposes the functions
under the `native` module. Mainly for testing purposes.

<!-- cargo-rdme end -->
