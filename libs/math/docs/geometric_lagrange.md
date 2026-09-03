# Exact geometric domains and interpolation {#math_geometric_lagrange}

@tableofcontents

`geometric_sequence_domain<FieldType, ValueType>` represents the exact-size domain

    x_i = r^i,  0 <= i < m,

where `r` is the field's configured geometric generator. The domain provides a linear-time field-element API for
evaluating every Lagrange basis polynomial at an arbitrary point and a context-based API for exact geometric
interpolation.

The relevant headers are:

| Facility | Header |
|---|---|
| Montgomery batch inversion | `<nil/crypto3/math/algorithms/batch_inverse.hpp>` |
| Evaluation-domain interface | `<nil/crypto3/math/domains/evaluation_domain.hpp>` |
| Exact geometric domain | `<nil/crypto3/math/domains/geometric_sequence_domain.hpp>` |
| Polynomial context | `<nil/crypto3/math/polynomial/backends/polynomial_backend.hpp>` |
| Schoolbook backend | `<nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>` |
| Mixed-radix backend | `<nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>` |

## Domain precomputation

Construction requires `m > 1`, a nonzero geometric generator, and distinct points `1, r, ..., r^(m-1)`. It rejects a
size that would repeat a point.

For the Lagrange path, the constructor precomputes:

* every domain point `x_i`;
* the barycentric weights `w_i = 1 / Z'(x_i)`; and
* the coefficients of the vanishing polynomial

      Z(X) = product(X - x_i, i = 0 .. m - 1).

For exact geometric interpolation, it also precomputes:

* the triangular powers `T_i = r^(i * (i - 1) / 2)` and their inverses;
* the denominator products `D_0 = 1` and `D_i = product(1 - r^j, j = 1 .. i)`;
* the inverses `1 / D_i`; and
* the fixed interpolation kernel `T_i / D_i`.

For general points, constructing every `w_i` from pairwise differences would take quadratic work. Geometric points
instead satisfy the recurrence

    w_0 = product((1 - x_j)^-1, j = 1 .. m - 1),
    w_i = -w_(i-1) * x_(m-1-i)^-1 * (1 - x_(m-i)) / (1 - x_i).

All required inverses of `1 - r^i`, together with `r^-1`, are obtained through one Montgomery batch inversion. For
`m` nonzero inputs, `batch_inverse_nonzero` performs one field inversion and exactly `3 * (m - 1)` field
multiplications.

Writing

    Z(X) = sum(c_j * X^(m-j), j = 0 .. m),  c_0 = 1,

the vanishing-polynomial coefficients use the recurrence

    c_j = -c_(j-1) * r^(j-1) * (1 - r^(m-j+1)) / (1 - r^j).

If `r` has exact order `m`, the final denominator vanishes and the implementation uses the resulting special case
`Z(X) = X^m - 1`.

Construction therefore takes `O(m)` field operations, one field inversion, and `O(m)` stored field elements.

## Exact geometric interpolation

The concrete geometric domain provides a non-virtual, compile-time backend-selected API:

```cpp
template<polynomial_arithmetic::PolynomialBackend Backend>
typename Backend::polynomial_type interpolate(
    const std::vector<typename Backend::polynomial_type::value_type>& evaluations,
    polynomial_arithmetic::polynomial_context<Backend>& context
) const;
```

The input must contain exactly one evaluation for each point `1, r, ..., r^(m-1)`. Any other count throws
`std::invalid_argument`. The domain points remain `FieldType::value_type`, while the evaluations and returned
coefficients may belong to a compatible extension field such as Fq12.

Interpolation proceeds in seven stages:

1. Validate that the input contains exactly `m` evaluations.
2. Build `scaled[i] = evaluations[i] * (-1)^i / D_i` and embed the fixed kernel `T_i / D_i` into the backend's
   coefficient field.
3. Multiply `scaled` by the fixed kernel through the supplied context and retain coefficients `0` through `m-1`.
4. Recover each Newton coefficient by multiplying convolution coefficient `i` by `1 / T_i`, then form the dynamic
   Newton input by multiplying it by `D_i`.
5. Form the fixed Newton-to-monomial kernel `(-1)^i * T_i / D_i` and embed it in reverse order.
6. Multiply the reversed fixed kernel by the dynamic Newton input through the same supplied context.
7. Set output coefficient `i` to product coefficient `m - 1 + i` multiplied by `1 / D_i`, then remove trailing
   zeros.

The second product is transposed multiplication: the fixed Newton-to-monomial kernel is reversed, not the dynamic
Newton input.

A schoolbook context requires no transform configuration:

```cpp
using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
polynomial_arithmetic::polynomial_context<backend_type> context;

const auto coefficients = domain.interpolate(evaluations, context);
```

A mixed-radix context must use a valid transform order supporting at least `2 * m - 1` product coefficients:

```cpp
using backend_type = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
polynomial_arithmetic::polynomial_context<backend_type> context {backend_type(transform_order)};

const auto coefficients = domain.interpolate(evaluations, context);
```

Both products use the caller's context. Interpolation does not construct a backend, invoke the legacy transform, or
perform field inversions.

## Evaluating all weights

The combined overload returns both the weights and `Z(t)` without repeating the product:

```cpp
#include <nil/crypto3/math/domains/geometric_sequence_domain.hpp>

namespace math = nil::crypto3::math;

math::geometric_sequence_domain<field_type> domain(domain_size);

value_type vanishing_at_t;
const std::vector<value_type> weights =
    domain.evaluate_all_lagrange_polynomials(t, vanishing_at_t);
```

Away from the domain, it uses the barycentric formula

    L_i(t) = Z(t) * w_i / (t - x_i).

The implementation forms all denominators `t - x_i`, computes `Z(t)` as their product, and batch-inverts the
denominators with one field inversion. If `t = x_i`, it instead returns the corresponding unit vector and sets
`Z(t) = 0` without attempting an inversion.

For evaluations `f(x_i)`, the returned weights satisfy

    f(t) = sum(f(x_i) * L_i(t), i = 0 .. m - 1)

for every polynomial of degree below `m`. The weights belong to `FieldType::value_type`; they may also scale values in
a compatible extension field.

## Rough current costs

| Operation | Rough current cost |
|---|---|
| Construct a size-`m` domain | `O(m)` field operations and one inversion |
| Interpolate `m` evaluations | Two backend products, `O(m)` additional field operations, and no inversions |
| Evaluate all weights at one off-domain point | `O(m)` field operations and one inversion |
| Evaluate at a domain point | `O(m)` work and no inversion |
| Evaluate all weights at `k` independent points | `O(k * m)` field operations and at most `k` inversions |

Returning `m` weights already requires linear output work. All scratch storage used by an evaluation or interpolation
is local, while the domain precomputation is immutable. Concurrent interpolation calls may share one domain, but each
call requires a separate polynomial context because backends may reuse mutable plans, caches, or scratch storage.
