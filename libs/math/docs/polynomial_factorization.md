# Polynomial factorization {#math_polynomial_factorization}

@tableofcontents

Crypto3.Math factors univariate coefficient polynomials over finite fields. The algorithms operate on the polynomial
type supplied by a multiplication backend, so the same factorization code can use the reference schoolbook backend or
a faster backend suitable for the coefficient field and operand sizes.

## Result representation

A complete factorization has the form

    input = leading_coefficient * product(factor.polynomial ^ factor.multiplicity).

`polynomial_factorization_result<Polynomial>` stores the original leading coefficient, the monic irreducible factors
and their positive multiplicities, and a `complete` flag. Nonconstant factors are canonical and monic. The canonical
zero polynomial and nonzero constants produce no polynomial factors; their scalar value is returned as
`leading_coefficient`.

Distinct-degree factorization uses a related result type. Each `distinct_degree_factor` contains the product of all
irreducible factors having one degree. Such a group is generally not itself irreducible. For example, a group reported
with `irreducible_factor_degree == 2` may be the product of several distinct irreducible quadratics.

## Factorization pipeline

Complete factorization composes three stages:

| Stage | Public API | Purpose |
|---|---|---|
| Square-free factorization | `square_free_factorization` | Separate factors by their multiplicity using Yun's algorithm. |
| Distinct-degree factorization | `distinct_degree_factorization_kaltofen_shoup` | Group the square-free factors by irreducible degree using blocked Frobenius steps and GCDs. |
| Equal-degree factorization | `equal_degree_factorization` | Split one degree group into individual irreducible factors using Cantor-Zassenhaus. |

`distinct_degree_factorization_reference` provides the unblocked distinct-degree algorithm as a correctness reference.
Production callers normally use the Kaltofen-Shoup implementation selected by `complete_factorization`.

The [polynomial arithmetic infrastructure](@ref math_polynomial_arithmetic) provides the reusable divisor contexts,
Newton division, GCD and half-GCD, quotient-ring operations, Brent-Kung modular composition, and iterated Frobenius
maps used by these stages.

### Header map

| Facility | Header |
|---|---|
| Result and callback types | `<nil/crypto3/math/polynomial/polynomial_factorization.hpp>` |
| Square-free factorization | `<nil/crypto3/math/polynomial/square_free_factorization.hpp>` |
| Reference distinct-degree factorization | `<nil/crypto3/math/polynomial/distinct_degree_factorization.hpp>` |
| Kaltofen-Shoup distinct-degree factorization | `<nil/crypto3/math/polynomial/kaltofen_shoup_distinct_degree_factorization.hpp>` |
| Equal-degree factorization | `<nil/crypto3/math/polynomial/equal_degree_factorization.hpp>` |
| Complete factorization | `<nil/crypto3/math/polynomial/complete_factorization.hpp>` |

## Complete factorization

The simplest entry point takes a polynomial arithmetic context and a caller-owned random generator:

```cpp
#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/math/polynomial/complete_factorization.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>
#include <nil/crypto3/random/algebraic_engine.hpp>

namespace math = nil::crypto3::math;
namespace pa = math::polynomial_arithmetic;
namespace fields = nil::crypto3::algebra::fields;

using field_type = fields::babybear;
using value_type = field_type::value_type;
using backend_type = pa::schoolbook_backend<value_type>;
using polynomial_type = backend_type::polynomial_type;

// (X + 1) * (X + 2)^2 = X^3 + 5 X^2 + 8 X + 4.
polynomial_type input = {
    value_type(4), value_type(8), value_type(5), value_type::one()
};

pa::polynomial_context<backend_type> arithmetic_context;
nil::crypto3::random::algebraic_engine<field_type> generator(42);
auto result = math::complete_factorization<backend_type>(
    input, arithmetic_context, generator);

for (const auto &factor : result.factors) {
    // factor.polynomial is monic and irreducible.
    // factor.multiplicity is its multiplicity in input.
}
```

The generator is injected rather than created internally. It must return independent, uniformly distributed elements
of the coefficient field. Supplying a seeded generator makes a run reproducible.

Cantor-Zassenhaus repeats randomized split trials until a nontrivial factor is found; the implementation does not set
an arbitrary retry limit. A generator that does not adequately explore the coefficient field can therefore prevent
progress. Callers that require a time or work limit should enforce it through their generator or at a higher level.

For large products, callers can instead construct a `mixed_radix_backend<RootField, ValueType>` with an adequate
transform order and place it in the polynomial context. This permits base-field roots to transform extension-field
coefficients. Backend choice changes the arithmetic implementation, not the factorization API or result.

## Staged factorization

The factorization APIs have callback overloads for callers that can decide after each emitted factor whether more work
is necessary:

```cpp
auto result = math::complete_factorization<backend_type>(
    input, arithmetic_context, generator,
    [](const math::polynomial_factor<polynomial_type> &factor) {
        return should_stop(factor)
            ? math::factorization_control::stop_factorization
            : math::factorization_control::continue_factorization;
    });
```

The factor that triggers the stop is included in `result.factors`, and `result.complete` is `false`. A complete run has
`complete == true`. Square-free and distinct-degree factorization also expose staged callbacks, allowing a caller to
stop at the earliest useful stage rather than computing a complete irreducible factorization.

## Preconditions and current restrictions

* Inputs are nonempty coefficient polynomials. Algorithms canonicalize trailing zero coefficients and return monic
  nonconstant factors.
* Square-free factorization currently requires the coefficient-field characteristic to exceed the input degree. This
  avoids the separate polynomial p-th-root path required when differentiation erases p-th powers.
* Distinct-degree and equal-degree entry points require square-free input. Complete factorization establishes this
  precondition through its first stage.
* Equal-degree and complete factorization currently implement the odd-characteristic Cantor-Zassenhaus algorithm and
  reject characteristic two.
* An equal-degree input must be a group whose irreducible factors all have the supplied degree. The API checks whether
  the factor degree divides the total degree but does not repeat distinct-degree factorization to verify the full
  precondition.

## Verifying a result

A complete result can be checked by starting with the constant polynomial containing `leading_coefficient`, multiplying
each monic factor into it `multiplicity` times through the same backend, and comparing the canonical result with the
input. This reconstruction identity is the primary contract of `polynomial_factorization_result`.

## Cost overview

Let `n` be the input degree, `d` the common irreducible-factor degree in an equal-degree group, and `Q` the
coefficient-field order. The current implementations perform the rough work shown below.

The table counts expensive polynomial operations. Their individual bounds—particularly GCD, exact division, modular
multiplication, and cached Frobenius maps—are given in the
[polynomial arithmetic cost overview](@ref math_polynomial_arithmetic) and depend on the selected multiplication
backend.

| Stage | Rough current work |
|---|---|
| Yun square-free factorization | Up to `n` iterations, each with one GCD and two exact divisions; divisor inverses are rebuilt |
| Reference distinct-degree factorization | Up to `n / 2` cached Frobenius maps and `n / 2` GCDs |
| Kaltofen-Shoup coarse phase | `O(sqrt(n))` Frobenius maps, `O(n)` modular multiplications, and `O(sqrt(n))` coarse GCDs |
| Kaltofen-Shoup fine splitting | Up to `O(n)` additional GCDs and exact divisions |
| One Cantor-Zassenhaus trial | `O(d * log(Q))` modular multiplications plus at most two polynomial GCDs |
| Equal-degree factorization | An expected linear number of Cantor-Zassenhaus trials in the number of output factors |
| Complete factorization | The sum of the square-free, Kaltofen-Shoup, and equal-degree stages |

The Cantor-Zassenhaus bounds are expected bounds for uniform random samples; a caller-supplied generator has no fixed
retry limit. Each Frobenius map uses the current Brent-Kung implementation and therefore includes its quadratic
coefficient-combination phase. Staged callbacks stop paying these costs once the caller has obtained enough factors.
