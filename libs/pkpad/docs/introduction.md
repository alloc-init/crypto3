# Introduction # {#pkpad_introduction}

@tableofcontents

Crypto3.Pkpad provides encoding and verification primitives used by public-key
schemes. The current headers include EME encodings such as OAEP and PKCS #1,
EMSA encodings, MGF1, and ISO/IEC 9796 support.

Crypto3.Pkpad consists of several parts to review:
* [Manual](@ref pkpad_manual).
* [Implementation](@ref pkpad_impl).
* [Concepts](@ref pkpad_concepts).

## Dependencies ## {#pkpad_dependencies}

Internal dependencies:

1. Crypto3.Codec
2. Crypto3.Hash
3. Crypto3.Algebra
4. Crypto3 algebra marshalling

Outer dependencies:
1. [Boost](https://boost.org)
