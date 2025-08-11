# Introduction # {#modes_introduction}

@tableofcontents

The Crypto3.Modes library extends the =nil; Foundation's cryptography suite and provides a set of block and
 stream cipher modes implemented in way C++ standard library implies: concepts, algorithms, predictable
  behavior, latest standard features support and clean architecture without compromising security and performance.
 
Crypto3.Modes consists of several parts to review:
* [Manual](@ref modes_manual).
* [Implementation](@ref modes_impl).
* [Concepts](@ref modes_concepts).

## Dependencies ## {#modes_dependencies}

Internal dependencies:

1. [Crypto3.Block](https://github.com/alloc-init/block.git)
2. [Crypto3.Stream](https://github.com/alloc-init/stream.git)
3. [Crypto3.Codec](https://github.com/alloc-init/codec.git)

Outer dependencies:
1. [Boost](https://boost.org) (>= 1.58)