# Introduction # {#codec_introduction}

@tableofcontents

The Crypto3.Codec library extends the =nil; Crypto3 C++
 cryptography suite and provides a set of encoding
algorithms implemented in way C++ standard library
implies: concepts, algorithms, predictable behavior,
latest standard features support and clean architecture
without compromising security and performance.
  
Crypto3.Codec consists of several parts to review:
* [Manual](@ref codec_manual).
* [Implementation](@ref codec_impl).
* [Concepts](@ref codec_concepts).

## Dependencies ## {#codec_dependencies}

In-suite dependencies:

Base58 support conditionally depends on Crypto3.Multiprecision. Other codecs do
not link another Crypto3 component directly.

Outer dependencies:
1. [Boost](https://boost.org)
