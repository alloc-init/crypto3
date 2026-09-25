# Marshalling utilities for [[alloc] init]'s Zero-Knowledge Schemes Cryptography

This module extends the [Crypto3 marshalling](../) utilities for
[Crypto3 zero-knowledge schemes](../../zk/).

## Building

This library is built as part of the Crypto3 monorepo. Follow the
[repository build instructions](../../../README.md#clone-and-build); the shared
[CMake modules](../../../cmake/modules/) are the repository's only submodule.

## Modified SAP SNARK

The [modified SAP example](example/modified_sap.cpp) converts an ordinary R1CS
circuit, generates proving and verification keys using Crypto3 ChaCha as the
setup randomness source, and serializes and deserializes the constraint system,
keys and proof. It proves one assignment with the restored proving key, then
verifies the restored proof against the caller's expected public input using
the restored verification key.

From the repository root:

```sh
cmake -S . -B build -DBUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target marshalling_modified_sap_example --parallel
./build/libs/marshalling/zk/example/marshalling_modified_sap_example
```

The example's fixed seed is for reproducibility; its setup keys are insecure.
See the [serialization specification](../../zk/docs/modified_sap_snark.md#9-serialization)
for the `fill_*` / `write` / `read` / `make_*` APIs, supported big-endian BN254
profile, integer-width compatibility, validation and input-size requirements.

## Dependencies

### Internal

* [Crypto3.Multiprecision](../../multiprecision/)
* [Crypto3.Algebra](../../algebra/)
* [Crypto3.ZK](../../zk/)
* [Crypto3.Marshalling Core](../core/)
* [Crypto3.Multiprecision Marshalling](../multiprecision/)
* [Crypto3.Algebra Marshalling](../algebra/)

### External

* [Boost](https://boost.org) (>= 1.74)
