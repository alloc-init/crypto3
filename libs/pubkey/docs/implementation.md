# Implementation {#pubkey_impl}

@tableofcontents

## Architecture {#pubkey_arch}

Crypto3.Pubkey separates scheme policies, key objects, processing modes,
accumulators, and generic algorithms. Scheme policies define cryptographic
types and operations; key wrappers retain key material and expose the operation
modes supported by a scheme.

The current public processing modes are isomorphic operation and verifiable
encryption. Secret-sharing schemes, including Shamir, Feldman, Pedersen, and
weighted Shamir, are provided independently under `pubkey/secret_sharing`.

## Algorithms {#pubkey_algorithms}

Public algorithms are declared under `nil/crypto3/pubkey/algorithm`. Available
operations include signing, verification, encryption, decryption, signature
aggregation, aggregate verification, and secret-sharing operations where the
selected scheme supports them.

Range overloads provide the shortest interface:

```cpp
signature_type signature = sign(message, private_key);
bool valid = verify(message, signature, public_key);
```

Iterator and output-iterator overloads are available for callers that need
explicit range or destination control.

## Keys And Policies {#pubkey_objects}

`private_key<Scheme>` and `public_key<Scheme>` adapt a scheme's native key
types to the generic algorithms. The private-key wrapper provides access to its
corresponding public key where the scheme supports deriving it.

Scheme definitions such as `bls`, `ecdsa`, and `eddsa` bind curve, hash, and
scheme-version policies into the concrete key and signature types.

## Accumulators {#pubkey_accumulators}

Algorithms feed input through a processing mode backed by a Boost.Accumulators
set. The mode selects the operation policy and result type. Range-returning
overloads own this state internally; overloads accepting an accumulator allow
incremental input before result extraction.
