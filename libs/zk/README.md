# Crypto3.ZK

Zero-knowledge proof systems, polynomial commitment schemes, arithmetization,
and supporting primitives for the Crypto3 suite.

Crypto3.ZK is a header-only C++ library. Its public API is under
[`include/nil/crypto3/zk`](include/nil/crypto3/zk/) in the
`nil::crypto3::zk` namespace, with proof-system and arithmetization types in
`nil::crypto3::zk::snark`. Link the CMake target `crypto3::zk`.

## Library Overview

| Area | Contents |
| --- | --- |
| [Algorithms](include/nil/crypto3/zk/algorithms/) | Generic `generate`, `prove`, `verify`, and `aggregate` adapters for proof systems. |
| [Polynomial commitments](include/nil/crypto3/zk/commitments/polynomial/) | FRI, LPC, KZG (including batched, v2, and IPP2 variants), and Pedersen commitments; Powers of Tau, proof-of-knowledge, proof-of-work, and MPC support. |
| [R1CS proof systems](include/nil/crypto3/zk/snark/systems/ppzksnark/) | Pairing-based R1CS ppzkSNARK and Groth16-style R1CS GG-ppzkSNARK implementations. The GG implementation includes basic, proof-aggregation, and encrypted-input modes. |
| [Placeholder](include/nil/crypto3/zk/snark/systems/plonk/placeholder/) | PLONK Placeholder preprocessing, proving, verification, permutation, gate, lookup, and dFRI argument components, with LPC/FRI and KZG commitment backends. |
| [Arithmetization](include/nil/crypto3/zk/snark/arithmetization/) | R1CS constraints; QAP and SAP arithmetic programs; PLONK variables, constraints, gates, lookup tables, assignment tables, and polynomial tables. |
| [Reductions](include/nil/crypto3/zk/snark/reductions/) | R1CS-to-QAP and R1CS-to-SAP instance and witness mappings. |
| [Transcripts](include/nil/crypto3/zk/transcript/) | Sequential and accumulative Fiat-Shamir transcripts for byte- and field-oriented hash functions. |
| [Expressions](include/nil/crypto3/zk/math/) | Constraint expressions, visitors, DAG conversion, cached evaluation, and permutation utilities. |
| [Routing](include/nil/crypto3/zk/snark/routing/) | Benes and AS-Waksman routing-network construction and validation. |

The library does not provide a single umbrella header. Include the header for
the scheme or component being used. The R1CS proof systems have scheme-level
headers, for example:

```cpp
#include <nil/crypto3/zk/snark/systems/ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <nil/crypto3/zk/algorithms/generate.hpp>
#include <nil/crypto3/zk/algorithms/prove.hpp>
#include <nil/crypto3/zk/algorithms/verify.hpp>
```

Given a compatible curve and an R1CS constraint system, the generic proof API
has the following shape:

```cpp
using proof_system = nil::crypto3::zk::snark::r1cs_gg_ppzksnark<curve_type>;

auto keypair = nil::crypto3::zk::generate<proof_system>(constraint_system);
auto proof = nil::crypto3::zk::prove<proof_system>(
    keypair.first, primary_input, auxiliary_input);
bool valid = nil::crypto3::zk::verify<proof_system>(
    keypair.second, primary_input, proof);
```

Complete demonstrations are available for [FRI](example/fri.cpp),
[KZG](example/kzg.cpp), and [Fiat-Shamir transcripts](example/transcript.cpp).
The [tests](test/) contain end-to-end examples for the R1CS and Placeholder
proof systems, commitment schemes, reductions, expressions, and routing.

## Building

Crypto3.ZK requires CMake 3.22 or newer and a C++23-capable compiler. Configure
the Crypto3 monorepo and build the ZK test target with:

```sh
cmake -S . -B build -DBUILD_TESTS=ON
cmake --build build --target zk_tests --parallel
ctest --test-dir build --output-on-failure
```

See the [root build instructions](../../README.md#clone-and-build) for platform
requirements and dependency setup. To consume the library from a CMake project
that includes Crypto3:

```cmake
add_subdirectory(path/to/crypto3)
target_link_libraries(my_target PRIVATE crypto3::zk)
```

The interface target links the Crypto3 Containers, Math, Algebra, Hash, and
Multiprecision components together with Boost.Container and Boost.Log. ZK tests
also use Crypto3.Random, Crypto3 marshalling, and Boost.Test.
