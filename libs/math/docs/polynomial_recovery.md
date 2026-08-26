# Polynomial recovery {#math_polynomial_recovery}

@tableofcontents

Crypto3.Math provides two complementary recovery primitives:

* square testing and square roots in a finite polynomial quotient field; and
* bounded rational reconstruction from a residue modulo a polynomial.

Together they can recover representations of irreducible factors by the polynomial norm form

    P(X)^2 - X * Q(X)^2.

The library exposes the generic arithmetic operations. Applications remain responsible for selecting factors, degree
bounds, random sources, and any policy for combining or rejecting recovered representations.

### Header map

| Facility | Header |
|---|---|
| Field orders and multiplicative-group decomposition | `<nil/crypto3/algebra/fields/field_order.hpp>` |
| Coefficient-field square-root helpers | `<nil/crypto3/algebra/fields/field_algorithms.hpp>` |
| Quotient-field square testing and roots | `<nil/crypto3/math/polynomial/polynomial_square_root.hpp>` |
| Bounded rational reconstruction | `<nil/crypto3/math/polynomial/polynomial_rational_reconstruction.hpp>` |

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

1. `polynomial_divisor_context<Backend>` caches the divisor and the inverse needed for modular reduction.
2. `polynomial_square_root_context<Backend>` caches the decomposition
   `order(K)^d - 1 = odd_order * 2^two_adicity` and a suitable quadratic nonresidue raised to `odd_order`.

```cpp
#include <algorithm>
#include <nil/crypto3/math/polynomial/polynomial_square_root.hpp>

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
#include <nil/crypto3/math/polynomial/polynomial_rational_reconstruction.hpp>

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

## Recovering a polynomial norm representation

The square-root and reconstruction APIs fit together as follows. For one irreducible factor `B` of degree `d`:

1. Represent the indeterminate by `x = {0, 1}` and test whether it is a square modulo `B`.
2. If it is square, compute `R` such that `R^2 = X mod B`.
3. Rationally reconstruct `R` with

       maximum_numerator_degree   = floor(d / 2),
       maximum_denominator_degree = floor((d - 1) / 2).

   This produces `P = R * Q mod B`.
4. Squaring the congruence gives

       P^2 - X * Q^2 = 0 mod B.

   The degree bounds make the left side have degree at most `d`, so it is a scalar multiple of `B`.

This is the local representation of `B` by the norm from adjoining a square root of `X`. Representations compose
multiplicatively:

    (P1^2 - X Q1^2) * (P2^2 - X Q2^2)
      = (P1 P2 + X Q1 Q2)^2 - X * (P1 Q2 + P2 Q1)^2.

Consequently, a caller can factor a target polynomial, recover eligible irreducible factors independently, account for
multiplicities and the scalar leading coefficient, and combine the local representations. Crypto3.Math deliberately
keeps that application-level policy separate from the generic factorization, quotient-field square-root, and bounded
rational-reconstruction primitives.

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

This example performs the same local steps required for a higher-degree factor: factor, test whether `X` is square in
the quotient, compute its root, reconstruct bounded `P` and `Q`, and normalize the remaining coefficient-field
scalar. Combining several eligible factors then uses the multiplicative identity above.

## Reuse and performance

The [polynomial arithmetic infrastructure](@ref math_polynomial_arithmetic) describes these contexts in detail. Reuse
one polynomial arithmetic context throughout a recovery, one divisor context for every operation modulo the same `B`,
and one square-root context for repeated roots modulo that divisor. This avoids rebuilding divisor inverses,
multiplicative-group decompositions, nonresidue powers, backend plans, and scratch storage.

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
