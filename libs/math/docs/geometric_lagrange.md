# Exact geometric-domain Lagrange weights {#math_geometric_lagrange}

@tableofcontents

`geometric_sequence_domain<FieldType, ValueType>` represents the exact-size domain

    x_i = r^i,  0 <= i < m,

where `r` is the field's configured geometric generator. The domain provides a linear-time field-element API for
evaluating every Lagrange basis polynomial at an arbitrary point.

The relevant headers are:

| Facility | Header |
|---|---|
| Montgomery batch inversion | `<nil/crypto3/math/algorithms/batch_inverse.hpp>` |
| Evaluation-domain interface | `<nil/crypto3/math/domains/evaluation_domain.hpp>` |
| Exact geometric domain | `<nil/crypto3/math/domains/geometric_sequence_domain.hpp>` |

## Domain precomputation

Construction requires `m > 1`, a nonzero geometric generator, and distinct points `1, r, ..., r^(m-1)`. It rejects a
size that would repeat a point.

For the Lagrange path, the constructor precomputes:

* every domain point `x_i`;
* the barycentric weights `w_i = 1 / Z'(x_i)`; and
* the coefficients of the vanishing polynomial

      Z(X) = product(X - x_i, i = 0 .. m - 1).

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
| Evaluate all weights at one off-domain point | `O(m)` field operations and one inversion |
| Evaluate at a domain point | `O(m)` work and no inversion |
| Evaluate all weights at `k` independent points | `O(k * m)` field operations and at most `k` inversions |

Returning `m` weights already requires linear output work. All scratch storage used by an evaluation is local, while
the domain precomputation is immutable, so concurrent weight evaluations on one domain do not share mutable scratch.
