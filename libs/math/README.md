# Crypto3.Math

Matrix and polynomial arithmetic, evaluation domains, and Fast Fourier
Transform algorithms for the Crypto3 suite.

This header-only component is maintained as part of the Crypto3 monorepo. See
the [root build instructions](../../README.md#clone-and-build) to configure the
project and run its tests. Link the component as `crypto3::math` from CMake.

API notes are available under [`docs`](docs/), and buildable examples are under
[`example`](example/).

## Matrix API

The header-only matrix API is under
`include/nil/crypto3/math/matrix`. Applications should use the Crypto3
frontends rather than access a backend directly:

- `regular_matrix<T>` and `regular_vector<T>` provide dense storage;
- `compressed_matrix<T>` and `compressed_vector<T>` provide sparse storage;
- matrix and vector expressions remain lazy until materialized into a frontend;
- `product`, `inner_product`, `element_product`, and row/column/subrange views
  are available from `operators.hpp`, while the frontends provide compound
  arithmetic;
- `find_element` and `for_each_nonzero` provide sparse access without exposing
  backend iterators;
- `solve.hpp` provides backend-generic Gaussian elimination,
  `back_substitute`, and matrix-shape helpers.

The currently supplied storage aliases use Boost.uBLAS, which is header-only.
`crypto3::math` publishes `Boost::headers` as an interface dependency, so
consumers only need to link `crypto3::math`:

```cmake
target_link_libraries(my_target PRIVATE crypto3::math)
```

`solve` accepts owning Crypto3 `matrix<Backend>` and `vector<Backend>` values,
uses the supplied backends as its working storage, and returns the same vector
backend. Non-owning expression views should be materialized into the desired
frontend before calling it.
