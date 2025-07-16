# Introduction # {#passhash_introduction}

@tableofcontents

The Crypto3.Passhash library extends the =nil; Foundation's cryptography suite and provides a set password hashes
 implemented in way C++ standard library implies: concepts, algorithms, predictable behavior, latest standard features support and clean architecture without compromising security and performance.
 
Crypto3.Passhash consists of several parts to review:
* [Manual](@ref passhash_manual).
* [Implementation](@ref passhash_impl).
* [Concepts](@ref passhash_concepts).

## Dependencies ## {#passhash_dependencies}

Internal dependencies:

1. [Crypto3.Mac](https://github.com/alloc-init/block.git)
2. [Crypto3.Codec](https://github.com/alloc-init/codec.git)
3. [Crypto3.Pbkdf](https://github.com/alloc-init/pbkdf.git)

Outer dependencies:
1. [Boost](https://boost.org) (>= 1.58)