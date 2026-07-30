# Manual {#zk_manual}

@tableofcontents

## Integration

Link the header-only library through CMake:

```cmake
target_link_libraries(my_target PRIVATE crypto3::zk)
```

There is no library-wide umbrella header. Include the concrete proof system,
commitment scheme, algorithms, and transcript types used by the application.

## R1CS GG-ppzkSNARK

The generic adapters provide the setup, proving, and verification flow:

```cpp
#include <nil/crypto3/zk/algorithms/generate.hpp>
#include <nil/crypto3/zk/algorithms/prove.hpp>
#include <nil/crypto3/zk/algorithms/verify.hpp>
#include <nil/crypto3/zk/snark/systems/ppzksnark/r1cs_gg_ppzksnark.hpp>

using proof_system =
    nil::crypto3::zk::snark::r1cs_gg_ppzksnark<curve_type>;

auto keypair = nil::crypto3::zk::generate<proof_system>(constraint_system);
auto proof = nil::crypto3::zk::prove<proof_system>(
    keypair.first, primary_input, auxiliary_input);
bool valid = nil::crypto3::zk::verify<proof_system>(
    keypair.second, primary_input, proof);
```

`constraint_system` is a
`snark::r1cs_constraint_system<curve_type::scalar_field_type>`. Set its
`primary_input_size` and `auxiliary_input_size`, then add
`snark::r1cs_constraint` objects. Each constraint contains three linear
combinations. Variable index `0` is constant one, indices
`1..primary_input_size` are public, and the remaining indices are private.

Check the inputs before setup:

```cpp
if (!constraint_system.is_valid() ||
    !constraint_system.is_satisfied(primary_input, auxiliary_input)) {
    throw std::invalid_argument("invalid or unsatisfied R1CS instance");
}
```

Applications that verify many proofs with one key can process the verification
key once and pass the resulting processed key to the same `verify` adapter. See
[`run_r1cs_gg_ppzksnark.hpp`](../test/systems/ppzksnark/r1cs_gg_ppzksnark/run_r1cs_gg_ppzksnark.hpp)
for both verification paths.

## Placeholder

Placeholder accepts a PLONK constraint system, assignment table, table
description, and commitment scheme. After defining `field_type` and
`placeholder_params_type`, the protocol sequence is:

```cpp
using public_preprocessor =
    nil::crypto3::zk::snark::placeholder_public_preprocessor<
        field_type, placeholder_params_type>;
using private_preprocessor =
    nil::crypto3::zk::snark::placeholder_private_preprocessor<
        field_type, placeholder_params_type>;
using prover = nil::crypto3::zk::snark::placeholder_prover<
    field_type, placeholder_params_type>;
using verifier = nil::crypto3::zk::snark::placeholder_verifier<
    field_type, placeholder_params_type>;

auto public_data = public_preprocessor::process(
    constraint_system, assignments.public_table(), table_description,
    prover_scheme);
auto private_data = private_preprocessor::process(
    constraint_system, assignments.private_table(), table_description);

auto proof = prover::process(
    public_data, std::move(private_data), table_description,
    constraint_system, prover_scheme);

commitment_scheme_type verifier_scheme(commitment_params);
bool valid = verifier::process(
    *public_data.common_data, proof, table_description, constraint_system,
    verifier_scheme);
```

Use a fresh verifier commitment-scheme instance configured with the same
parameters; do not reuse mutable prover state. The complete LPC/FRI and KZG
configurations are maintained in
[`placeholder_test_runner.hpp`](../test/systems/plonk/placeholder/placeholder_test_runner.hpp).

## Polynomial commitments

FRI uses the algorithm adapters directly:

```cpp
auto tree = nil::crypto3::zk::algorithms::precommit<fri_type>(
    polynomial, params.D[0], params.step_list[0]);
auto commitment = nil::crypto3::zk::algorithms::commit<fri_type>(tree);
auto proof = nil::crypto3::zk::algorithms::proof_eval<fri_type>(
    polynomial, tree, params, prover_transcript);
bool valid = nil::crypto3::zk::algorithms::verify_eval<fri_type>(
    proof, commitment, params, verifier_transcript);
```

The prover and verifier transcripts must be initialized identically. See the
[FRI example](../example/fri.cpp) for coefficient and DFS polynomial forms.

The [KZG example](../example/kzg.cpp) demonstrates single-polynomial and
batched commitments. Its setup is generated from a known secret for testing;
production applications must use securely generated parameters.

LPC is stateful. Add polynomials to batches, commit, register evaluation
points, and produce a proof. Configure a fresh verifier with matching batch
sizes and evaluation points before calling `verify_eval()`. See
[`test/commitment/lpc.cpp`](../test/commitment/lpc.cpp) for the complete order.

## Transcripts

The recommended transcript is
`transcript::fiat_shamir_heuristic_sequential<HashType>`:

```cpp
#include <nil/crypto3/zk/transcript/fiat_shamir.hpp>

std::vector<std::uint8_t> initialization = {0, 1, 2, 3};
nil::crypto3::zk::transcript::fiat_shamir_heuristic_sequential<hash_type>
    transcript(initialization);

transcript(value_to_absorb);
auto challenge = transcript.challenge<field_type>();
```

The verifier must repeat initialization, absorption, and challenge requests in
the same order. The [transcript example](../example/transcript.cpp) shows
Keccak and Poseidon configurations.

## Examples and tests

The `example/` directory contains maintained source demonstrations for FRI,
KZG, and transcripts. These sources are not currently added to the standard
top-level CMake build. End-to-end, buildable usage is available in `test/`.

Build and run the ZK tests from the repository root:

```sh
cmake -S . -B build -DBUILD_TESTS=ON
cmake --build build --target zk_tests --parallel
ctest --test-dir build --output-on-failure
```
