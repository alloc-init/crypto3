# Manual # {#fft_manual}

@tableofcontents

Crypto3.Math provides coefficient-form and evaluation-form polynomials, evaluation domains, FFT algorithms, and
backend-aware polynomial arithmetic. The higher-level polynomial algorithms share a
`polynomial_arithmetic::polynomial_context<Backend>`. Reusing a context lets its multiplication backend retain plans,
configuration, and scratch storage across an operation.

The following pages describe the optimized domain and polynomial facilities:

* [Exact geometric-domain Lagrange weights](@ref math_geometric_lagrange) covers exact-size geometric points,
  barycentric precomputation, batch inversion, vanishing polynomials, and linear-time weight evaluation.
* [Polynomial arithmetic infrastructure](@ref math_polynomial_arithmetic) describes multiplication backends,
  reusable contexts, division, GCD, quotient-ring arithmetic, modular composition, and Frobenius maps shared by the
  higher-level algorithms.
* [Polynomial factorization](@ref math_polynomial_factorization) covers square-free, distinct-degree, equal-degree,
  complete, and staged factorization.
* [Polynomial recovery](@ref math_polynomial_recovery) covers square testing and square roots in polynomial quotient
  fields, bounded rational reconstruction, and the relation between these operations and polynomial norms.
