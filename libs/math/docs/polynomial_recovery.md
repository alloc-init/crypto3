# Polynomial recovery {#math_polynomial_recovery}

@tableofcontents

Crypto3.Math provides three complementary recovery facilities:

* square testing and square roots in a finite polynomial quotient field;
* bounded rational reconstruction from a residue modulo a polynomial; and
* one-call reconstruction of the fixed `X`-norm representation of an irreducible polynomial.

Together they can recover representations of irreducible factors by the polynomial norm form

    P(X)^2 - X * Q(X)^2.

The one-call operation owns the degree bounds and exact normalization for this norm equation. The caller supplies the
polynomial arithmetic context and the coefficient-field generator.

### Header map

| Facility | Header |
|---|---|
| Field orders and multiplicative-group decomposition | `<nil/crypto3/algebra/fields/field_order.hpp>` |
| Coefficient-field square-root helpers | `<nil/crypto3/algebra/fields/field_algorithms.hpp>` |
| Quotient-field square testing and roots | `<nil/crypto3/math/polynomial/quotient_ring/polynomial_square_root.hpp>` |
| Bounded rational reconstruction | `<nil/crypto3/math/polynomial/reconstruction/polynomial_rational_reconstruction.hpp>` |
| One-call `X`-norm reconstruction | `<nil/crypto3/math/polynomial/reconstruction/polynomial_x_norm_reconstruction.hpp>` |

## Field orders and coefficient square roots

The quotient-field algorithms derive their exponents from the public field-order utilities:

* `field_characteristic<Field>()` returns the prime characteristic;
* `field_order<Field>()` returns the number of elements, including the extension degree described by `Field`;
* `extension_field_order<Field>(d)` returns `field_order<Field>()^d`; and
* the corresponding multiplicative-group decomposition functions write an order minus one as
  `odd_order * 2^two_adicity`.

The templates accept either a Crypto3 field type or its value type and return `boost::multiprecision::cpp_int`, since
an extension-field order can exceed the fixed-width field representation. These integers are used as public loop
bounds and exponents; polynomial coefficients remain native field values.

BN254 Fq12 values provide `is_square()` and `sqrt()`. `is_square()` handles zero and applies the finite-field square
criterion. `sqrt()` requires a square input and asserts that precondition; callers handling arbitrary values should
test first. The implementation uses the field's multiplicative-group decomposition and Tonelli-Shanks rather than
converting the value to a generic integer representation.

## Quotient-field square testing

Let `K` be a finite field and let `B` be irreducible of degree `d`. The quotient `K[X]/(B)` is a field with
`order(K)^d` elements. A polynomial with degree below `d` is the canonical representative of one quotient-field
element.

`is_square_mod(input, divisor_context, arithmetic_context)` reports whether that element is a square. It treats zero
as a square and uses Euler's criterion for a general nonzero representative. For the canonical indeterminate `X`, the
implementation can use

    Norm(X) = (-1)^d * B(0) / leading_coefficient(B)

when the coefficient value type provides `is_square()`. This replaces an expensive quotient-ring exponentiation with
one coefficient-field square test.

Irreducibility of `B` is a caller precondition. The function does not factor `B`, and Euler's criterion does not
characterize squares if the quotient has zero divisors.

## Quotient-field square roots

`square_root_mod` uses Tonelli-Shanks in `K[X]/(B)`. Repeated calls should share two immutable precomputations:

* `polynomial_divisor_context<Backend>` caches the divisor and the inverse needed for modular reduction.
* `polynomial_square_root_context<Backend>` caches the decomposition
  `order(K)^d - 1 = odd_order * 2^two_adicity` and a suitable quadratic nonresidue raised to `odd_order`.

```cpp
#include <algorithm>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_square_root.hpp>

// B must be irreducible. Its degree is d.
const std::size_t d = B.size() - 1;
const std::size_t inverse_precision = std::max<std::size_t>(1, d - 1);

pa::polynomial_context<backend_type> arithmetic_context;
math::polynomial_divisor_context<backend_type> divisor_context(
    B, inverse_precision, arithmetic_context);

// known_nonresidue is canonical, reduced modulo B, and nonsquare in K[X]/(B).
math::polynomial_square_root_context<backend_type> square_root_context(
    known_nonresidue, divisor_context, arithmetic_context);

polynomial_type root;
if (math::square_root_mod(root, value, square_root_context,
                          arithmetic_context)) {
    // root^2 == value mod B.
}
```

Instead of supplying a known nonresidue, a second context constructor accepts a caller-owned generator. The generator
returns canonical reduced polynomial representatives until the context finds a nonsquare. The divisor context must
outlive the square-root context that references it.

The search has no fixed retry limit. The generator must be capable of producing a quotient-field nonsquare; a
pathological generator can make context construction fail to progress. Applications that need a work limit should
enforce one in the supplied generator.

The square-root operation supports odd characteristic. It returns `false` and stores zero when the input is not a
square; zero is returned as its own root. Output may alias input. The coefficient field may itself be an extension
field, including BN254 Fq12.

## Bounded rational reconstruction

Given a reduced residue `R` modulo a nonconstant polynomial `B`, `rational_reconstruct` searches for polynomials `P`
and `Q` such that

    P = R * Q mod B,

subject to caller-supplied degree bounds. It follows the Euclidean remainder sequence until the remainder reaches the
numerator bound, while tracking the corresponding coefficient of `R`. The returned denominator is monic, and the
numerator is scaled by the same field element.

```cpp
#include <nil/crypto3/math/polynomial/reconstruction/polynomial_rational_reconstruction.hpp>

polynomial_type numerator;
polynomial_type denominator;
const bool recovered = math::rational_reconstruct(
    numerator, denominator, residue, modulus,
    maximum_numerator_degree, maximum_denominator_degree,
    arithmetic_context);
```

The strict condition

    maximum_numerator_degree + maximum_denominator_degree < degree(modulus)

ensures uniqueness. The function returns `false` if the Euclidean candidate exceeds the denominator bound and leaves
both outputs unchanged. It rejects noncanonical inputs, a constant modulus, an unreduced residue, nonunique bounds, or
using the same object for both outputs. Either output may otherwise alias an input.

For example, over any field in which the displayed small integers are distinct, the residue

    R = 2 + X + X^2 + X^3

modulo

    B = 3 + 3X + 3X^2 + 2X^3 + X^4

reconstructs with bounds `degree(P) <= 0` and `degree(Q) <= 2` as

    P = 1,
    Q = 2 + 2X + X^2.

## One-call `X`-norm reconstruction

`recover_polynomial_x_norm_representation` composes quotient-field square roots, bounded rational reconstruction, and
coefficient-field normalization to recover polynomials `P` and `Q` satisfying

    P^2 - X * Q^2 = g.

The name explicitly identifies `X` as the quadratic element: the represented element is `P + Q * sqrt(X)`, whose norm
is `P^2 - X * Q^2`. The input `g` must be a canonical nonconstant irreducible polynomial. Irreducibility is a caller
precondition and is not tested.

```cpp
#include <nil/crypto3/math/polynomial/reconstruction/polynomial_x_norm_reconstruction.hpp>

// arithmetic_context and coefficient_generator are caller-owned.
auto representation = math::recover_polynomial_x_norm_representation<backend_type>(
    g, arithmetic_context, coefficient_generator);

if (representation) {
    polynomial_type exact_norm = math::evaluate_polynomial_x_norm<backend_type>(
        *representation, arithmetic_context);
    // exact_norm == g
}
```

The coefficient generator returns coefficient-field values. The recovery operation adapts those values into canonical
degree-below-`degree(g)` representatives of `K[X]/(g)`. The generator remains caller-owned and must eventually supply
coefficients forming a quotient-field nonsquare.

On success, the optional contains `polynomial_x_norm_representation<polynomial_type>`. Its `p` and `q` members hold the
two recovered coefficient polynomials.

The function reduces `{0, 1}` modulo `g`, including when `g` is linear, and recovers `R` with `R^2 = X mod g`. It then
reconstructs `P = R * Q mod g` with

    degree(P) <= floor(degree(g) / 2),
    degree(Q) <= floor((degree(g) - 1) / 2).

These bounds imply

    P^2 - X * Q^2 = lambda * g

for a coefficient-field scalar `lambda`. A nonzero square `lambda` is removed by scaling both outputs by
`sqrt(lambda^-1)`. The normalized norm is evaluated again and compared exactly with `g` before success is returned.
Both polynomial squares use the caller's arithmetic context; multiplication by `X` is a coefficient shift.

| Outcome | Contract |
|---|---|
| Representation returned | The degree bounds hold and `evaluate_polynomial_x_norm(result, context) == g`. |
| No value returned | `X` is nonsquare modulo `g`, bounded reconstruction fails, or `lambda` is zero or nonsquare. |
| `std::invalid_argument` | The input is empty, noncanonical, zero, or constant, or a composed API contract is violated. |
| `std::logic_error` | An operation reported success but a required modular, scalar-multiple, or final exact identity is inconsistent. |

This is the local representation of `g` by the norm from adjoining a square root of `X`. Representations compose
multiplicatively:

    (P1^2 - X Q1^2) * (P2^2 - X Q2^2)
      = (P1 P2 + X Q1 Q2)^2 - X * (P1 Q2 + P2 Q1)^2.

Consequently, a caller can factor a target polynomial, recover eligible irreducible factors independently, account for
multiplicities and the scalar leading coefficient, and combine the local representations. Crypto3.Math deliberately
keeps that application-level policy separate from the generic factorization, quotient-field square-root, bounded
rational-reconstruction, and one-call `X`-norm reconstruction facilities.

### Worked example over F7

Consider

    H = X^2 + 1

over the field with seven elements. Its only possible roots would square to `-1 = 6`, but the squares in F7 are
`0, 1, 2, 4`, so `H` is irreducible and its complete factorization contains just `H`.

In the quotient field `F7[X]/(H)`, the polynomial

    R = 2 + 2X

is a square root of `X`:

    R^2 = 4 + 8X + 4X^2 = X mod H.

For `d = 2`, reconstruction uses numerator degree at most one and denominator degree at most zero. It returns

    P = 2 + 2X,
    Q = 1,

so that

    P^2 - X Q^2 = 4 + 7X + 4X^2 = 4H.

The reconstruction determines a representation up to this nonzero scalar. Since `4^-1 = 2` and `3^2 = 2` in F7,
scale both recovered polynomials by `3`:

    P' = 6 + 6X,
    Q' = 3.

Then the scalar is corrected exactly:

    P'^2 - X Q'^2 = H.

The one-call API performs this square test, quotient-field square root, bounded reconstruction, scalar normalization,
and final exact verification using the supplied polynomial context and coefficient generator. Combining several
eligible factors then uses the multiplicative identity above.

## Reuse and performance

The [polynomial arithmetic infrastructure](@ref math_polynomial_arithmetic) describes these contexts in detail. The
one-call API reuses the caller's polynomial arithmetic context throughout recovery. It constructs one divisor context
and one square-root context and reuses them for all operations within that call. Direct users of the lower-level APIs
can retain those contexts across repeated operations modulo the same divisor.

Let `Q` be the coefficient-field order, `d` the irreducible divisor degree, and `s` the two-adicity of `Q^d - 1`.
The current implementations have the rough bounds shown below.

The table counts quotient-ring multiplications and Euclidean steps. The cost of each underlying polynomial
multiplication, reduction, and division is given in the
[polynomial arithmetic cost overview](@ref math_polynomial_arithmetic) and depends on the selected backend.

| Operation | Rough current bound |
|---|---|
| General `is_square_mod` | `O(d * log Q)` quotient-ring multiplications |
| `is_square_mod(X)` norm shortcut | One coefficient-field square test |
| Square-root context construction | Expected `O(d * log Q)` quotient-ring multiplications |
| `square_root_mod` | `O(d * log Q + s^2)` quotient-ring multiplications |
| Rational reconstruction | Up to `O(d)` sequential Euclidean steps; `O(d^2)` with schoolbook short-by-long products |

The square-root context has no deterministic bound with a caller-supplied generator because it searches until it
finds a nonresidue. Rational reconstruction does not currently use half-GCD acceleration.
