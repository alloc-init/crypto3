# Polynomial arithmetic infrastructure {#math_polynomial_arithmetic}

@tableofcontents

The higher-level polynomial algorithms in Crypto3.Math share a backend-aware arithmetic layer. It separates the
mathematical algorithms from the multiplication implementation and gives callers explicit ownership of expensive
precomputation.

The main dependency chain is:

| Layer | Main types and operations |
|---|---|
| Multiplication | `PolynomialBackend`, `schoolbook_backend`, `mixed_radix_backend`, `bn254_small_prime_fft_backend` |
| Arithmetic lifetime | `polynomial_context` and `polynomial_context_options` |
| Division | `inverse_series`, `polynomial_divisor_context`, `divrem`, `remainder`, `exact_division` |
| GCD | `gcd`, with optional internal half-GCD reduction |
| Quotient ring | `mulmod`, `squaremod`, `powmod`, `powmod_x` |
| Modular composition | `compose_mod_reference`, `polynomial_composition_precomputation`, `compose_mod` |
| Frobenius | `polynomial_frobenius_context`, `frobenius_map`, `frobenius_product`, `frobenius_weighted_sum` |

### Header map

| Facility | Header |
|---|---|
| Representation concepts | `<nil/crypto3/math/polynomial/concepts.hpp>` |
| Backend concept and arithmetic context | `<nil/crypto3/math/polynomial/backends/polynomial_backend.hpp>` |
| Schoolbook and mixed-radix backends | `<nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>`, `<nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>` |
| BN254 Fp12 small-prime backend | `<nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>` |
| Power-series inversion | `<nil/crypto3/math/polynomial/arithmetic/power_series.hpp>` |
| Divisor context and division | `<nil/crypto3/math/polynomial/arithmetic/polynomial_division.hpp>` |
| GCD | `<nil/crypto3/math/polynomial/arithmetic/gcd.hpp>` |
| Quotient-ring multiplication | `<nil/crypto3/math/polynomial/quotient_ring/polynomial_modular_arithmetic.hpp>` |
| Quotient-ring exponentiation | `<nil/crypto3/math/polynomial/quotient_ring/polynomial_exponentiation.hpp>` |
| Modular composition | `<nil/crypto3/math/polynomial/quotient_ring/polynomial_composition.hpp>` |
| Frobenius maps, products, and weighted sums | `<nil/crypto3/math/polynomial/quotient_ring/polynomial_frobenius.hpp>` |

Internal headers live beside their owning subsystem. Relative to `nil/crypto3/math/`, the small-prime backend's
parameters, arithmetic, transforms, conversion, CRT export, and pointwise kernels are grouped under
`polynomial/backends/detail/small_prime_fft/`. The prepared composition helper lives in
`polynomial/quotient_ring/detail/`. General helpers shared across math components remain in `detail/`.

## Multiplication backends

`PolynomialBackend` is the compile-time interface used by the arithmetic layer. A backend supplies an associated
coefficient polynomial type and three alias-safe operations:

```cpp
backend.multiply(output, left, right);
backend.square(output, input);
backend.multiply_low(output, left, right, coefficient_count);
```

Inputs and outputs use ascending coefficient order and canonical representation. `multiply_low` computes the product
modulo `X^coefficient_count`; a zero coefficient count produces the canonical zero polynomial.

Crypto3.Math currently provides three implementations:

* `schoolbook_backend<ValueType>` uses direct quadratic coefficient multiplication. It is the reference backend and is
  usually preferable for small operands or products with one very short operand.
* `mixed_radix_backend<RootFieldType, ValueType>` uses roots from `RootFieldType` to transform `ValueType`
  coefficients. Its constructor takes an explicit maximum transform order, caches plans for that order's divisors,
  and reuses transform scratch storage. The configured order must cover every requested product.
* `bn254_small_prime_fft_backend` multiplies BN254 Fp12 coefficient polynomials through auxiliary small-prime
  transforms and exact signed CRT reconstruction. It takes a power-of-two maximum transform length and caches all
  smaller power-of-two plans. Scalar and linear operands use direct coefficient arithmetic.

This separation permits, for example, roots in a base field to transform values in an extension field. Backend
selection is explicit; higher-level algorithms do not silently replace the caller's backend.

### BN254 Fp12 small-prime products

This backend has a fixed coefficient type, available as `bn254_small_prime_fft_backend::value_type`. It converts
the twelve base-field coordinates of each coefficient into residues modulo auxiliary 50-bit primes. Each prime has
twelve transform channels; pointwise operations apply the BN254 tower identities in that prime's quotient ring.
Inverse transforms and signed CRT reconstruction return exact BN254 Fp12 coefficients.

For base-field modulus `p` and at most `t` overlapping polynomial terms, the conservative coordinate bound is

    C = 102 * t * (p - 1)^2.

The product of the selected primes must exceed `2*C`. Construction selects the shortest sufficient catalogue prefix
for the configured maximum length; all smaller plans share that prefix, coefficient converter, and CRT data. The
prime count is therefore fixed for a backend instance rather than selected again for every call. The 32,768-point
configuration uses eleven primes. CRT export uses a certified direct reduction into Montgomery representation when
possible and an exact reconstruction fallback otherwise.

The FMA kernels use integer-valued binary64 residues with proved bounds. They require the supported floating-point
compiler settings and round-to-nearest mode; other supported settings or rounding modes select exact integer
arithmetic. The calling thread's rounding mode must remain unchanged throughout an operation. AVX2/FMA and AVX-512
transform kernels depend on the compiler target, and IFMA acceleration also depends on the project's feature
configuration. A binary must run on a CPU supporting its compiled instruction set. This backend executes serially
and provides no constant-time guarantee.

The constructor is

```cpp
bn254_small_prime_fft_backend(max_transform_size,
                             max_spectrum_bytes,
                             max_matrix_spectrum_bytes);
```

Only `max_transform_size` is required. Ordinary FFT products need a transform containing the full linear product of
the retained input prefixes. `max_spectrum_bytes` defaults to the maximum `size_t` and bounds **one** spectrum at
the configured maximum length. A spectrum contains `prime_count * 12 * transform_length` doubles. Another working
spectrum, twiddle tables, output coefficients, and independently prepared operands require additional memory.
`max_matrix_spectrum_bytes` defaults to 256 MiB and bounds additional shared input spectra for optional matrix
products; zero disables that capability. Neither setting is a process-wide memory limit.

For example, an explicitly configured context can be constructed as follows:

```cpp
#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>

namespace pa = nil::crypto3::math::polynomial_arithmetic;
using backend_type = pa::bn254_small_prime_fft_backend;

pa::polynomial_context<backend_type> arithmetic_context(backend_type(32768));
```

Construction checks transform support, CRT capacity, and the spectrum limit before allocating plans. Product calls
reuse spectrum capacity, and dedicated squaring converts and transforms only one operand. Results are assembled
before replacing the output, preserving aliasing support and leaving output unchanged on failure.

## Polynomial contexts

`polynomial_arithmetic::polynomial_context<Backend>` owns one backend and the algorithm-selection parameters used by
division, GCD, composition, factorization, and recovery. Keeping the context alive lets a stateful backend reuse plans,
configuration, and scratch storage. Continuing with the small-prime backend above:

```cpp
namespace pa = nil::crypto3::math::polynomial_arithmetic;

pa::polynomial_context_options options;
options.gcd_half_gcd_cutoff = 0;

pa::polynomial_context<backend_type> arithmetic_context(
    backend_type(32768), options);
```

The context forwards `multiply`, `square`, and `multiply_low` to its backend. It is intended for sequential reuse.
Separate contexts are required for concurrent calls when the backend owns mutable scratch storage.

### Algorithm-selection options

| Option | Default | Meaning |
|---|---:|---|
| `basecase_divisor_coefficient_cutoff` | 10 | Use quadratic long division when the divisor has at most this many coefficients. |
| `basecase_quotient_coefficient_cutoff` | 2 | Use quadratic long division when the quotient has at most this many coefficients. |
| `half_gcd_basecase_cutoff` | 30 | At or below this coefficient count, recursive half-GCD constructs its transformation iteratively; zero disables this base case. |
| `gcd_half_gcd_cutoff` | 0 | Use half-GCD when the smaller operand reaches this size; zero disables half-GCD. |
| `modular_composition_cached_power_limit` | maximum `size_t` | Cap the number of Brent-Kung baby powers; constructing a composition cache requires a positive limit. |
| `distinct_degree_screening_limit` | 0 | Screen initial irreducible-factor degrees before the full Kaltofen-Shoup baby table; zero disables screening. |
| `modular_composition_prepared_byte_limit` | 0 | Bound the additional prepared Fp12 coefficient payload per composition cache; zero disables preparation. |
| `modular_composition_simd_lanes` | 8 | Select four or eight IFMA lanes when supported; zero selects scalar prepared dots. The byte budget still controls whether preparation is enabled. |

The two division cutoffs are independent: satisfying either enabled cutoff selects long division. Setting one of them
to zero disables that criterion. Half-GCD is disabled by default because its crossover depends strongly on the
coefficient type and multiplication backend; callers should choose a nonzero cutoff from representative benchmarks.

## Newton power-series inversion

`inverse_series(output, input, coefficient_count, arithmetic_context)` computes

    output * input = 1 mod X^coefficient_count.

The constant coefficient of `input` must be nonzero. Starting from its scalar inverse, Newton iteration doubles the
known precision on each step. The implementation preserves the already correct low block and computes only the new
high coefficients using `multiply_low`. Output may alias input, and a zero requested precision returns the canonical
zero polynomial.

Power-series inversion is the precomputation behind fast polynomial division.

## Polynomial divisor contexts

`polynomial_divisor_context<Backend>` stores a canonical nonzero divisor `B` and a truncated inverse of its reversed
polynomial. If `d = degree(B)`, define

    reverse(B) = X^d * B(X^-1).

Construction precomputes

    reverse(B)^-1 mod X^inverse_precision.

The chosen precision is the maximum number of quotient coefficients that later fast divisions may require. Reusing the
context avoids repeating this Newton inversion for every division or modular reduction by the same `B`.

```cpp
const std::size_t quotient_coefficient_bound =
    dividend.size() >= divisor.size()
        ? dividend.size() - divisor.size() + 1
        : 1;
math::polynomial_divisor_context<backend_type> divisor_context(
    divisor, quotient_coefficient_bound, arithmetic_context);
```

The context is immutable after construction. Operations receiving it assume it represents the intended divisor; when
another precomputation also depends on `B`, callers must keep those contexts paired with the same polynomial.

### Fixed-operand FFT preparation

With either FFT backend, construction can also prepare the FFT of the reversed-inverse prefix and one representation
of the divisor. These spectra own their data and belong to the divisor snapshot; changing the source polynomial
does not change them. The cache is bounded to at most two spectra and never grows during division. Constant divisors
and configurations selecting long division need no spectra.

The backend exposes optional `prepare_low_product` / `try_multiply_low_prepared` operations through the arithmetic
context. A preparation records the fixed prefix, precision and transform length. Reuse requires the same precision
and a compatible transform; the small-prime backend also checks the catalogue prefix. Otherwise the divisor context
uses ordinary `multiply_low`. Changing quotient lengths does not build more cache entries. Backends such as schoolbook
need only implement the three required operations.

The second spectrum is either the original low-product divisor prefix or, when a strictly shorter supported
transform exists, the complete divisor folded for cyclic remainder reconstruction. Both are never retained together.
Preparation pays one forward FFT per channel of each nonzero fixed operand up front, plus spectrum storage and
allocation. For the small-prime backend, the channels are the selected primes times twelve coordinates. Repeated
compatible operations then avoid those forward FFTs. This trades setup time and memory for lower repeated cost.

## Polynomial division

`divrem` computes canonical `quotient` and `remainder` satisfying

    dividend = quotient * B + remainder,
    degree(remainder) < degree(B).

The quotient and remainder must be distinct output objects, but either may alias the dividend. For small divisors or
quotients, the arithmetic options select quadratic long division. Otherwise, for

    k = degree(dividend) - degree(B) + 1,

reversal turns division into the truncated product

    reverse(quotient) = reverse(dividend) * reverse(B)^-1 mod X^k.

The divisor context must have inverse precision of at least `k` on this fast path. Long division does not use the
precomputed inverse and therefore does not require that precision.

### Short cyclic remainder reconstruction

Once the exact quotient is known, the identity `A = Q*B + R` gives, for a supported FFT order `M >= degree(B)`:

    R = (A mod (X^M - 1)) - (Q*B mod (X^M - 1)).

The high coefficients wrap identically in both terms and cancel. Since `degree(R) < degree(B) <= M`, the difference
is the exact remainder. The implementation folds the dividend and the **complete** divisor, including its leading
coefficient. In particular, that coefficient wraps into cell zero when `M == degree(B)`.

Both FFT backends use this path for `degree(A) <= 2*degree(B) - 2`, the range from multiplying reduced residues.
Each selects the smallest supported order at least `degree(B)` and uses it only when shorter than the ordinary
low-product transform. Outside that range, or when the active backend cannot reuse the prepared order, division
uses its original reconstruction. Quotient computation and inverse-precision validation are unchanged.

`prepare_cyclic_remainder` and `try_cyclic_remainder` are optional backend/context operations. The latter requires
the exact quotient for the same divisor used in preparation; it returns `false` without changing output when the
degree range or transform is incompatible. Its output may alias the dividend or quotient. The separate preparation
type prevents confusing a cyclic convolution with a low product.

Generic `multiply_low(A, B, k)` continues to mean `A*B mod X^k`: an ordinary k-point cyclic convolution cannot
replace it. Both reconstruction methods use the existing FFT kernels and preserve canonical results.

For dense `squaremod` calls using inverse-based division and compatible nonzero preparations, squaring, quotient
calculation and remainder reconstruction each need one forward and one inverse FFT per channel. Cyclic reconstruction
shortens only the remainder transforms; squaring and quotient calculation retain their existing transform lengths. Setup
performs the fixed-operand forward transforms once. These operation counts do not directly predict execution time:
pointwise products, folding, allocation and other division work also contribute. For the small-prime backend,
coefficient conversion and CRT export also contribute to each product.

Two convenience operations share the same dispatch:

* `remainder` returns only `dividend mod B`.
* `exact_division` returns the quotient and throws `std::invalid_argument` if the remainder is nonzero.

## GCD and half-GCD

`gcd(output, left, right, arithmetic_context)` returns the canonical monic greatest common divisor. Output may alias an
input, `gcd(0, B)` is the monic form of `B`, and `gcd(0, 0)` is zero.

The Euclidean algorithm repeatedly replaces `(A, B)` by `(B, A mod B)`. Each step preserves the common divisors and
strictly lowers the second degree until the final nonzero remainder is the GCD up to a scalar.

For sufficiently large operands, half-GCD batches several Euclidean steps. One quotient step acts on the polynomial
pair through

    [ 0   1 ] [ A ]   [ B         ]
    [ 1  -q ] [ B ] = [ A - q * B ].

Half-GCD recursively derives a product of these two-by-two polynomial matrices from the high coefficient halves, then
applies it to the complete inputs. A reduction lowers the second polynomial to roughly half the original first size.
The public `gcd` operation selects this internal path only when `gcd_half_gcd_cutoff` is nonzero and the smaller operand
reaches the cutoff. At or below `half_gcd_basecase_cutoff`, it constructs the same transformation iteratively.

Short-quotient Euclidean steps use direct division without constructing a reciprocal context. Recursive half-GCD
reuses the already transformed high halves when reconstructing the full pair. Backends may also expose
`try_multiply_matrix<Columns>` for a row-major 2-by-2 polynomial matrix times a 2-by-`Columns` matrix. The small-prime
backend supports one or two columns, sharing input spectra and reconstructing each sum of products once when its
shape, CRT bound, and workspace budget permit. Outputs must be distinct and pointers nonnull; outputs may alias
inputs. A declined operation leaves all outputs unchanged, and half-GCD uses ordinary products and additions.

Half-GCD reduces the number of sequential Euclidean divisions but introduces polynomial-matrix products and temporary
storage. It is therefore not unconditionally faster; its cutoff should be tuned for the active backend and coefficient
field.

## Quotient-ring arithmetic

The modular operations work in `K[X]/(B)` for any nonzero `B`. Irreducibility is not required unless a caller needs the
quotient to be a field.

* `mulmod` computes `left * right mod B`.
* `squaremod` uses the backend's dedicated square operation and reduces the result.
* `powmod` performs binary exponentiation, reducing the base and every intermediate product.
* `powmod_x` computes `X^exponent mod B` with the same binary schedule, replacing multiplication by the base with
  a coefficient shift and reduction of its single possible leading term.

Outputs are canonical and may alias their inputs. If `B` has degree `d`, multiplying representatives already reduced
modulo `B` requires at most `d - 1` inverse coefficients. Reducing an initially unreduced operand may require more.
Modulo a nonzero constant, every result is the zero polynomial because the quotient is the zero ring.

Both power functions accept built-in integer exponents and compatible multiprecision integer types. Exponents must
be nonnegative; exponent zero returns the quotient-ring identity. For a positive `L`-bit exponent, `powmod_x` performs
`L - 1` modular squares and at most `L - 1` linear multiply-by-X passes. It handles nonmonic and linear divisors.

## Brent-Kung modular composition

Modular composition computes

    outer(inner(X)) mod B.

`compose_mod_reference` uses Horner's rule and one modular multiplication per nonleading coefficient of `outer`. It is
the simple correctness reference.

The faster `compose_mod` overload uses blocked Brent-Kung composition. For block size `k`, write

    outer(Y) = F0(Y) + F1(Y) * Y^k + F2(Y) * Y^(2k) + ... .

`polynomial_composition_precomputation` caches the baby powers

    1, inner, ..., inner^(k - 1) mod B

and the giant step `inner^k mod B`. Each `Fi(inner)` is a linear combination of the baby powers, and the block values
are combined by Horner's rule in the giant step. The default `k = ceil(sqrt(L))`, where `L` is the maximum outer
coefficient count, balances precomputation against giant-step multiplications.

Construct and reuse the precomputation when composing several outer polynomials with the same inner polynomial and
divisor. The one-off overload constructs it internally. `modular_composition_cached_power_limit` caps `k`, reducing
memory from cached powers at the cost of more giant-step multiplications. Zero is rejected by the Brent-Kung APIs;
it does not select Horner composition automatically. `compose_mod_reference` remains available explicitly.

### Prepared coefficient combinations

Supported Fp12 coefficient types can additionally prepare the linear forms used by the block dot products. Set
`modular_composition_prepared_byte_limit` before constructing the composition cache to enable this preparation.
Its default is zero. `prepared_coefficient_bytes()` reports the retained payload. The budget applies to each cache
separately and excludes the ordinary baby powers, giant step, temporary block outputs, and per-call scratch.

The scalar prepared kernel uses 54 base-field product coordinates from the Fp12 tower. BN254 builds with
`CRYPTO3_HAS_AVX512IFMA` can use a 23-coordinate evaluation/interpolation representation with four or eight lanes.
`modular_composition_simd_lanes` accepts 0, 4, or 8; zero keeps the scalar prepared kernel. Unsupported SIMD builds
retain the scalar prepared path when it fits the budget. Unsupported coefficient types or a budget too small for
any preparation retain the ordinary coefficient loop.

Prepared evaluation batches up to four adjacent blocks, prepares each selected outer coefficient once per block,
and accumulates base-field dot products before recombination. Blocks with fewer than eight nonzero nonconstant
terms use the ordinary loop. SIMD preparation additionally requires a cached coefficient width of at least four.
All paths retain descending Horner order and exact field results. Prepared data is owned by the precomputation and
read-only during composition; copies own independent data, and scratch belongs to the call. This optimization is
independent of the polynomial multiplication backend and can be used with schoolbook, mixed-radix, or small-prime products.

The formation of block values still costs `O(L * degree(B))` field operations and can dominate at large comparable
degrees. Brent-Kung reduces modular multiplications; it does not remove the quadratic work needed to form these linear
combinations. Prepared dots and SIMD improve its constants without changing that asymptotic bound.

## Iterated Frobenius maps

Let the finite coefficient field `K` contain `Q` elements. Every coefficient satisfies `a^Q = a`, so for a polynomial
`A`:

    A(X)^Q mod B = A(X^Q mod B) mod B.

`polynomial_frobenius_context<Backend>` computes `X^Q mod B` using `powmod_x`, owns the required divisor context, and
builds a Brent-Kung precomputation for composition with that value. The context can then apply many Frobenius maps
without repeating the field-order exponentiation or rebuilding baby powers.

```cpp
math::polynomial_frobenius_context<backend_type> frobenius_context(
    divisor, arithmetic_context);

math::frobenius_map(output, input, frobenius_context, arithmetic_context);
math::frobenius_map(output, input, iteration_count,
                    frobenius_context, arithmetic_context);
```

One map raises a quotient-ring element to the `Q`-th power. The iterated overload applies that map
`iteration_count` times, producing the `Q^iteration_count` power; zero iterations only reduce the input. The divisor
need not be irreducible. This is a polynomial quotient-ring Frobenius operation and is distinct from any specialized
coordinate-level Frobenius implementation supplied by an extension-field element type.

### Frobenius products and weighted sums

Two operations compute aggregates of conjugates by binary doubling. Let `F` be a power of the coefficient-field
Frobenius, so `F` fixes coefficients, and supply the canonical reduced image `F(X) = X^(Q^s) mod B`:

```cpp
math::frobenius_product(output, input, term_count, image,
                       divisor_context, arithmetic_context);
math::frobenius_weighted_sum(output, input, weight, term_count, image,
                            divisor_context, arithmetic_context);
```

The product is `input * F(input) * ... * F^(term_count - 1)(input)`. For the weighted sum, define

    P_0 = 1,
    P_i = weight * F(weight) * ... * F^(i - 1)(weight).

The result is `sum(P_i * F^i(input), i = 0 .. term_count - 1)`. An empty product is the quotient-ring identity;
an empty sum is zero. A singleton is the reduced input. Outputs are canonical and may alias any polynomial input.
Inputs must be canonical, and the supplied image must be reduced. Its algebraic identity is a caller precondition;
the functions do not repeat the exponentiation to verify it. Irreducibility of `B` is not required.

These aggregates use `O(log(term_count))` compositions and modular products, with composition-cache construction
included at each needed doubling. Only the initial cache and one doubling cache are live at once. The iterated
`frobenius_map` overload above still applies one map per iteration; it does not use this aggregate algorithm.

Square-root recovery also uses an internal characteristic-power product for extension coefficients. If `Q = p^e`,
raising to `p` requires a coefficient Frobenius twist as well as composition with `X^p`. The public aggregate APIs
fix the coefficient field; passing `X^p` to them in place of `X^Q` is generally incorrect when `e > 1`. See
[polynomial recovery](@ref math_polynomial_recovery) for the characteristic-sized square-root seed.

## Context ownership and reuse

The arithmetic context normally has the longest lifetime. Divisor, composition, Frobenius, and square-root contexts
are immutable precomputations tied to particular polynomials and may reference or logically depend on one another.
Construct them once per fixed input and reuse them sequentially. Do not combine precomputations built for different
divisors merely because their degrees match.

## Cost overview

Let `n` be the polynomial size, `M(n)` the active backend's multiplication cost, `e` an exponent, and `Q` the
coefficient-field order. These are rough bounds for the current implementation:

| Operation | Rough time bound |
|---|---|
| Schoolbook multiplication or squaring | `O(n^2)` |
| Mixed-radix multiplication | `O(N * sum(radices))`, where `N` is the selected transform size |
| Small-prime multiplication | `O(r * N * log N)` small-residue transform work, plus conversion, tower pointwise work, and CRT/export; `r` is the selected prime count and there are twelve coordinate channels per prime |
| Newton series inversion | `O(M(n))` |
| Divisor-context construction | `O(M(n))` |
| Long division | `O(n^2)` |
| Newton division with a reused divisor context | `O(M(n))` |
| Euclidean GCD | `O(n^2)` with the default small-quotient path |
| Half-GCD | `O(M(n) * log n)` |
| `mulmod` or `squaremod` | `O(M(n))` with fast division; `O(n^2)` with schoolbook arithmetic |
| `powmod` | `O(M(n) * log e)` with fast division |
| `powmod_x` | `O(M(n) * log e)` with fast division, with linear shifts/reductions for set bits |
| Reference modular composition | `O(n * M(n))` |
| Current Brent-Kung composition | `O(sqrt(n) * M(n) + n^2)` |
| Frobenius-context construction | `O(M(n) * log Q + sqrt(n) * M(n))` |
| One cached Frobenius map | `O(sqrt(n) * M(n) + n^2)` |
| Frobenius product or weighted sum of `t` terms | `O(log t)` compositions, composition-cache preparations, and modular products for `t > 1` |

The `n^2` term in Brent-Kung and Frobenius maps is the block-formation phase, including when prepared dots are enabled.
Both FFT backends' `multiply_low` operations truncate their inputs but still transform their complete prefix product,
so they have the same rough asymptotic cost as multiplying those prefixes in full. They do not currently use truncated
transforms. The transform-only cost for the small-prime backend is not a complete product-cost estimate: conversion
and exact export can be substantial and must be included in measurements.
