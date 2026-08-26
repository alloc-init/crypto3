# Polynomial arithmetic infrastructure {#math_polynomial_arithmetic}

@tableofcontents

The higher-level polynomial algorithms in Crypto3.Math share a backend-aware arithmetic layer. It separates the
mathematical algorithms from the multiplication implementation and gives callers explicit ownership of expensive
precomputation.

The main dependency chain is:

| Layer | Main types and operations |
|---|---|
| Multiplication | `PolynomialBackend`, `schoolbook_backend`, `mixed_radix_backend` |
| Arithmetic lifetime | `polynomial_context` and `polynomial_context_options` |
| Division | `inverse_series`, `polynomial_divisor_context`, `divrem`, `remainder`, `exact_division` |
| GCD | `gcd`, with optional internal half-GCD reduction |
| Quotient ring | `mulmod`, `squaremod`, `powmod` |
| Modular composition | `compose_mod_reference`, `polynomial_composition_precomputation`, `compose_mod` |
| Frobenius | `polynomial_frobenius_context`, `frobenius_map` |

### Header map

| Facility | Header |
|---|---|
| Representation concepts | `<nil/crypto3/math/polynomial/concepts.hpp>` |
| Backend concept and arithmetic context | `<nil/crypto3/math/polynomial/polynomial_backend.hpp>` |
| Schoolbook and mixed-radix backends | `<nil/crypto3/math/polynomial/schoolbook_backend.hpp>`, `<nil/crypto3/math/polynomial/mixed_radix_backend.hpp>` |
| Power-series inversion | `<nil/crypto3/math/polynomial/power_series.hpp>` |
| Divisor context and division | `<nil/crypto3/math/polynomial/polynomial_division.hpp>` |
| GCD | `<nil/crypto3/math/polynomial/gcd.hpp>` |
| Quotient-ring multiplication | `<nil/crypto3/math/polynomial/polynomial_modular_arithmetic.hpp>` |
| Quotient-ring exponentiation | `<nil/crypto3/math/polynomial/polynomial_exponentiation.hpp>` |
| Modular composition | `<nil/crypto3/math/polynomial/polynomial_composition.hpp>` |
| Frobenius maps | `<nil/crypto3/math/polynomial/polynomial_frobenius.hpp>` |

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

Crypto3.Math currently provides two implementations:

* `schoolbook_backend<ValueType>` uses direct quadratic coefficient multiplication. It is the reference backend and is
  usually preferable for small operands or products with one very short operand.
* `mixed_radix_backend<RootFieldType, ValueType>` uses roots from `RootFieldType` to transform `ValueType`
  coefficients. Its constructor takes an explicit maximum transform order, caches plans for that order's divisors,
  and reuses transform scratch storage. The configured order must cover every requested product.

This separation permits, for example, roots in a base field to transform values in an extension field. Backend
selection is explicit; higher-level algorithms do not silently replace the caller's backend.

## Polynomial contexts

`polynomial_arithmetic::polynomial_context<Backend>` owns one backend and the algorithm-selection parameters used by
division, GCD, composition, factorization, and recovery. Keeping the context alive lets a stateful backend reuse plans,
configuration, and scratch storage:

```cpp
namespace pa = nil::crypto3::math::polynomial_arithmetic;

pa::polynomial_context_options options;
options.gcd_half_gcd_cutoff = 0;

pa::polynomial_context<backend_type> arithmetic_context(
    backend_type{}, options);
```

The context forwards `multiply`, `square`, and `multiply_low` to its backend. It is intended for sequential reuse.
Separate contexts are required for concurrent calls when the backend owns mutable scratch storage.

### Algorithm-selection options

| Option | Default | Meaning |
|---|---:|---|
| `basecase_divisor_coefficient_cutoff` | 10 | Use quadratic long division when the divisor has at most this many coefficients. |
| `basecase_quotient_coefficient_cutoff` | 2 | Use quadratic long division when the quotient has at most this many coefficients. |
| `half_gcd_basecase_cutoff` | 30 | Below this size, recursive half-GCD constructs its transformation iteratively. |
| `gcd_half_gcd_cutoff` | 0 | Use half-GCD when the smaller operand reaches this size; zero disables half-GCD. |
| `modular_composition_cached_power_limit` | maximum `size_t` | Limit the number of Brent-Kung baby powers retained in memory. |

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
reaches the cutoff. Below `half_gcd_basecase_cutoff`, it constructs the same transformation iteratively.

Half-GCD reduces the number of sequential Euclidean divisions but introduces polynomial-matrix products and temporary
storage. It is therefore not unconditionally faster; its cutoff should be tuned for the active backend and coefficient
field.

## Quotient-ring arithmetic

The modular operations work in `K[X]/(B)` for any nonzero `B`. Irreducibility is not required unless a caller needs the
quotient to be a field.

* `mulmod` computes `left * right mod B`.
* `squaremod` uses the backend's dedicated square operation and reduces the result.
* `powmod` performs binary exponentiation, reducing the base and every intermediate product.

Outputs are canonical and may alias their inputs. If `B` has degree `d`, multiplying representatives already reduced
modulo `B` requires at most `d - 1` inverse coefficients. Reducing an initially unreduced operand may require more.
Modulo a nonzero constant, every result is the zero polynomial because the quotient is the zero ring.

`powmod` accepts built-in integer exponents and compatible multiprecision integer types. Exponents must be
nonnegative; exponent zero returns the quotient-ring identity.

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
memory from cached powers at the cost of more giant-step multiplications.

The current coefficientwise formation of block values costs `O(L * degree(B))` field operations and can dominate at
large comparable degrees. Brent-Kung reduces modular multiplications; it does not remove the quadratic work needed to
form these linear combinations.

## Iterated Frobenius maps

Let the finite coefficient field `K` contain `Q` elements. Every coefficient satisfies `a^Q = a`, so for a polynomial
`A`:

    A(X)^Q mod B = A(X^Q mod B) mod B.

`polynomial_frobenius_context<Backend>` computes and stores `X^Q mod B`, owns the required divisor context, and builds a
Brent-Kung precomputation for composition with that value. The context can then apply many Frobenius maps without
repeating the field-order exponentiation or rebuilding baby powers.

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
| Newton series inversion | `O(M(n))` |
| Divisor-context construction | `O(M(n))` |
| Long division | `O(n^2)` |
| Newton division with a reused divisor context | `O(M(n))` |
| Euclidean GCD | `O(n^2)` with the default small-quotient path |
| Half-GCD | `O(M(n) * log n)` |
| `mulmod` or `squaremod` | `O(M(n))` with fast division; `O(n^2)` with schoolbook arithmetic |
| `powmod` | `O(M(n) * log e)` with fast division |
| Reference modular composition | `O(n * M(n))` |
| Current Brent-Kung composition | `O(sqrt(n) * M(n) + n^2)` |
| Frobenius-context construction | `O(M(n) * log Q + sqrt(n) * M(n))` |
| One cached Frobenius map | `O(sqrt(n) * M(n) + n^2)` |

The `n^2` term in Brent-Kung and Frobenius maps is the current coefficient-by-coefficient block-formation phase.
Mixed-radix `multiply_low` truncates its inputs but still transforms their complete prefix product, so it has the same
rough asymptotic cost as multiplying those prefixes in full.
