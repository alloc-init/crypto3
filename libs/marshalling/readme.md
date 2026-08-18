# [[alloc] init]'s Marshalling library

This library is used throughout the project to transform data from one type to
another. To define representation rules for custom type we use template-defined
pseudo-DSL.

With term *marshalling* we denote both serialization and deserialization of
objects.

Serialization is the process of transforming C++ object (elliptic curve point
for example) into the series of simple units (vector of bytes for example).

Deserialization is the opposing process.

The process is done using intermediate representation of C++ objects as
tree-like structure of tuples/bundles/arrays of elements. Then this structure
is transformed into series of bytes.

For serialization the corresponding functions are `make_XXX` and `pack_XXX`.
Deserialization is done using functions like `fill_XXX` and `unpack_XXX`.

The library consists of following modules:

```
marshalling
├── algebra
├── containers
├── core
├── math
├── multiprecision
├── zk
```

For convenience the umbrella `pack` function is provided, that can be used for
either process. See `algebra/example` for example usage. `core/example`
contains samples for creating custom `fill_S` and `make_S` functions.

## Algebra

Provides support for marshalling of algebraic elements: field elements and
elliptic curve points.

## Containers

Provides support for marshalling of Merkle tree structures.

## Core

Provides support for marshalling of basic types as well as template-defined
helper functions for arrays and vectors of other elements.

## Math

Provides support for marshalling of math objects: polynomials, matrices, and 
arithmetic expressions.

## Multiprecision

Provides support for marshalling of multiprecision integers.

## ZK

Provides support for marshalling of Zero-Knowledge structures: commitments of
different schemes, assignment table, constraint system and various structures
comprising Placeholder proof.

## Build and Test

Configure the Crypto3 monorepo from its root, build the marshalling tests, and
run them through CTest:

```sh
cmake -S . -B build -DBUILD_TESTS=ON
cmake --build build --target tests --parallel
ctest --test-dir build --output-on-failure
```

To build and run only the algebra curve-element test:

```sh
cmake --build build --target marshalling_algebra_curve_element_test
ctest --test-dir build --output-on-failure \
  -R '^marshalling_algebra_curve_element_test$'
```

See the [root build instructions](../../README.md#clone-and-build) for
dependencies and platform-specific setup.

## Usage

The suite is used as header-only libraries.
