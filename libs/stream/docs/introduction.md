# Introduction # {#stream_ciphers_introduction}

@tableofcontents

The Crypto3.Stream library extends the Nil Foundation's cryptography suite and provides a set of stream ciphers
 implemented in way C++ standard library implies: concepts, algorithms, predictable behavior, latest standard features
  support and clean architecture without compromising security and performance.
  
Crypto3.Stream consists of several parts to review:
* [Manual](@ref stream_ciphers_manual).
* [Implementation](@ref stream_ciphers_impl).
* [Concepts](@ref stream_ciphers_concepts).
 
A small part of each topic is reviewed right at introduction, but the detailed information is recommended to look inside of a corresponding chapter.
   
## Algorithms ## {#stream_ciphers_algorithms}

Crypto3.Stream contains the following stream cipher:

* `stream::chacha`

## Dependencies ## {#stream_dependencies}

Internal dependencies:

1. Crypto3.Block

Outer dependencies:
1. [Boost](https://boost.org)
