# Introduction # {#functional_introduction}

@tableofcontents

Crypto3.Functional library extends [[alloc] init]'s Crypto3 C++ cryptography suite and provides a set of functional
encryption schemes implemented in a way C++
standard library implies: concepts, algorithms, predictable behavior, latest standard features support and clean
architecture without compromising security and performance.

Crypto3.Functional consists of several parts to review:

* [Manual](@ref functional_manual).
* [Implementation](@ref functional_impl).
* [Concepts](@ref functional_concepts).

## Dependencies ## {#functional_dependencies}

Internal dependencies:

1. [Crypto3.Block](https://github.com/alloc-init/block.git)
2. [Crypto3.Codec](https://github.com/alloc-init/codec.git)

External dependencies:

1. [Boost](https://boost.org) (>= 1.58)

## Contents

### Attribute based encryption

This package includes schemes for key policy (KP) and ciphertext policy (CP) ABE schemes. The KP-ABE scheme GPSW is
based on the paper by Goyal, Pandey, Sahai, Waters: "Attribute-Based Encryption for Fine-Grained Access
Control of Encrypted Data".

The CP-ABE scheme FAME is based on the paper by Shashank Agrawal and Melissa Chase: "FAME: Fast Attribute-based Message
Encryption".

### Functional encryption schemes for inner products.

Based on security assumptions, the schemes are organized into submodules simple (s-IND-CPA security), and fullysec (
e.g. "fully secure", offering adaptive IND-CPA security).

Note that in both modules you will find single input as well as multi input schemes. Construction of all multi input
schemes is based on the work of Abdalla et. al (see paper: https://eprint.iacr.org/2017/972.pdf)

Most of the functions write their results to pointers which are passed as parameters. The memory used by these results
needs to be freed by the caller with *xyz_free* functions (each struct defined in this library has a respective
*struct_name_free* function). These free the memory used by the struct's members, but not the memory used by the struct
itself. This is due to the initialization of the structs' members in functions such as *generate_keys* and
*derive_fe_key* - this makes the API simpler and easier to use as the user
of the library does not to know all struct fields.

The results of the functions are the parameters listed _before_ the pointer to the struct in the parameter list.
Consult the appropriate documentation for each scheme for more thorough descriptions of functions and their results.

#### Fully secure schemes for functional encryption of inner products.

All implementations in this package are based on the reference paper by Agrawal, Libert and Stehlé (
see https://eprint.iacr.org/2015/608.pdf), and offer adaptive security under chosen-plaintext attacks (IND-CPA
security).

The reference scheme is public key, which means that no master secret key is required for the encryption.

For instantiation from the decisional Diffie-Hellman assumption (DDH), see struct damgard (and its multi-input variant
damgard_multi, which is a secret key scheme, because a part of the secret key is required for the encryption).

For instantiation from learning with errors (LWE), see struct lwe_fs.

#### Simple schemes for functional encryption of inner products.

All implementations in this package are based on the reference paper by Abdalla et. al (see https://eprint.iacr.org/2015/017.pdf). The reference scheme offers selective security under chosen-plaintext attacks (s-IND-CPA security).

The reference scheme is public key, which means that no master secret key is required for the encryption.

For instantiation from the decisional Diffie-Hellman assumption (DDH), see struct ddh (and its multi-input variant ddh_multi, which is a secret key scheme, because a part of the secret key is required for the encryption).

For instantiation from learning with errors (LWE), see structs lwe and ring_lwe.