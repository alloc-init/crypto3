# C++ concepts {#fft_concepts}

@tableofcontents

This page lists the C++20 concepts exposed by Crypto3.Math. A concept is a compile-time predicate used to constrain a
template argument. It describes the interface an argument must provide; it does not construct an object or generate
runtime code.

Crypto3.Math currently exposes concepts for matrix and vector expressions, polynomial representations, and polynomial
multiplication backends. Concepts declared in a `detail` namespace are implementation constraints and are summarized
separately at the end of this page.

## Matrix and vector concepts

The matrix concepts are declared in:

    <nil/crypto3/math/matrix/concepts.hpp>

### Readable expressions

`VectorExpression<T>` describes a readable vector-like expression. The type must provide `value_type`, `size_type`,
`size()`, and indexed access through `value(i)`.

`MatrixExpression<T>` describes a readable matrix-like expression. The type must provide `value_type`, `size_type`,
indexed access through `value(row, column)`, and dimensions through either:

* `rows()` and `columns()`; or
* `size1()` and `size2()`.

The free `rows(expression)` and `columns(expression)` functions provide one interface for both dimension naming
conventions.

Expression concepts are deliberately weak. In particular, an expression-template result may be readable without
owning storage or permitting mutation.

### Writable backends

`VectorBackend<T>` extends `VectorExpression<T>`. It additionally requires a semiregular type and writable indexed
access through `value(i)`.

`MatrixBackend<T>` similarly extends `MatrixExpression<T>` with semiregular value semantics and writable access
through `value(row, column)`.

The resizable variants add their corresponding resize operation:

| Concept | Additional operation |
|---|---|
| `ResizableVectorBackend<T>` | `value.resize(size)` |
| `ResizableMatrixBackend<T>` | `value.resize(rows, columns)` |

Algorithms that only read their operands should accept expression types. Algorithms that create or overwrite results
use backend or resizable-backend constraints as appropriate.

## Polynomial representation concepts

The polynomial representation concepts are declared in:

    <nil/crypto3/math/polynomial/concepts.hpp>

Both concepts require `value_type`, `size_type`, `size()`, `degree()`, indexed read access, and a representation tag.
They accept const and reference-qualified types.

### CoefficientPolynomial

`CoefficientPolynomial<T>` identifies a polynomial stored as coefficients in ascending degree order:

    [a0, a1, ..., an] represents a0 + a1 X + ... + an X^n.

Its `representation_type` must be `coefficient_representation`. Crypto3's `polynomial`, `polynomial_view`, and
`polymorphic_polynomial` satisfy this concept.

```cpp
template<nil::crypto3::math::CoefficientPolynomial Polynomial>
void consume_coefficients(const Polynomial &polynomial);
```

### EvaluationPolynomial

`EvaluationPolynomial<T>` identifies a polynomial stored as evaluations over a domain. Indexed entries are samples,
not coefficients, and the logical degree may differ from `size() - 1`.

Its `representation_type` must be `evaluation_representation`. Crypto3's `polynomial_dfs`, `polynomial_dfs_view`, and
`polymorphic_polynomial_dfs` satisfy this concept.

```cpp
template<nil::crypto3::math::EvaluationPolynomial Polynomial>
void consume_evaluations(const Polynomial &polynomial);
```

The representation tag prevents a coefficient algorithm from accidentally accepting an evaluation polynomial merely
because both types provide similar container operations. A bare `std::vector` satisfies neither concept because it
does not identify what its entries represent.

These concepts describe readable representation only. They do not require field-valued entries, mutability,
canonical storage, arithmetic operators, or alias-safe operations.

## PolynomialBackend

The polynomial multiplication backend concept is declared in:

    <nil/crypto3/math/polynomial/backends/polynomial_backend.hpp>

`polynomial_arithmetic::PolynomialBackend<Backend>` requires an associated `polynomial_type` that satisfies
`CoefficientPolynomial` and three operations:

```cpp
backend.multiply(output, left, right);
backend.square(output, input);
backend.multiply_low(output, left, right, coefficient_count);
```

The operations produce canonical coefficient polynomials and permit output to alias an input. `multiply_low` computes
the product modulo `X^coefficient_count`.

Backend operations receive a mutable backend object so implementations can reuse plans and scratch storage. The
[polynomial arithmetic infrastructure](@ref math_polynomial_arithmetic) documents the supplied schoolbook and
mixed-radix backends and the `polynomial_context` that owns them.
