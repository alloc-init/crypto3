# Implementation {#zk_impl}

Crypto3.ZK is an `INTERFACE` library. Template implementations reside in the
public headers under `include/nil/crypto3/zk`; there is no library-wide umbrella
header.

@tableofcontents

## Public API layout

| Directory | Responsibility |
| --- | --- |
| `algorithms/` | Generic proof-system adapters: `generate`, `prove`, `verify`, and `aggregate`. |
| `commitments/` | Commitment traits, polynomial batching, FRI, LPC, KZG variants, Pedersen, and setup utilities. |
| `math/` | Linear combinations, expressions, visitors, DAG evaluation, and permutations. |
| `snark/arithmetization/` | R1CS, QAP, SAP, and PLONK data structures. |
| `snark/reductions/` | R1CS-to-QAP and R1CS-to-SAP mappings. |
| `snark/systems/` | R1CS ppzkSNARKs and PLONK Placeholder. |
| `snark/routing/` | Benes and AS-Waksman routing networks. |
| `transcript/` | Fiat-Shamir transcript implementations. |

## R1CS proof systems

`snark::r1cs_ppzksnark` and `snark::r1cs_gg_ppzksnark` are scheme facades. The
generic functions in `zk/algorithms` delegate to each facade's static
`generate`, `prove`, and `verify` functions. Keys and proofs are scheme-specific
types; the GG keypair is a `std::pair` whose `first` member is the proving key
and whose `second` member is the verification key.

The GG implementation separates generator, prover, strong/weak verifier, and
verification-key processing policies. `snark::proving_mode` selects basic,
aggregate, or encrypted-input specializations. The scheme-level header
`r1cs_gg_ppzksnark.hpp` includes those specializations.

R1CS-to-QAP and R1CS-to-SAP reductions expose `instance_map`,
`instance_map_with_evaluation`, and `witness_map` operations. The proof systems
use these mappings to convert constraints and assignments into polynomial
relations.

## Placeholder pipeline

Placeholder is assembled from cooperating templates rather than a single
scheme facade:

- `snark::placeholder_circuit_params` describes circuit column counts.
- `snark::placeholder_params` combines circuit and commitment parameters.
- `snark::placeholder_public_preprocessor` commits to fixed circuit data and
  creates verifier `common_data`.
- `snark::placeholder_private_preprocessor` converts witness columns.
- `snark::placeholder_prover` creates `snark::placeholder_proof`.
- `snark::placeholder_verifier` checks the proof, including permutation, gate,
  lookup, and commitment-opening arguments.

Transcript initialization binds the constraint system, table dimensions,
usable rows, commitment parameters, field parameters, and an application
identifier. The implementation is in
`snark/systems/plonk/placeholder/detail/transcript_initialization_context.hpp`.

## Commitment implementations

FRI exposes stateless algorithm adapters such as `precommit`, `commit`,
`proof_eval`, and `verify_eval`. LPC combines FRI with the stateful
`commitments::polys_evaluator` batching interface. KZG provides single and
batched variants, while `kzg_commitment_scheme` and
`kzg_commitment_scheme_v2` adapt batched KZG to Placeholder's expected
interface.

Commitment schemes intentionally use their own parameter, commitment, proof,
and transcript types. Include the concrete scheme header rather than depending
on transitive includes.

## Validation sources

The `test/` tree is the authoritative source for complete protocol assembly:

- `test/systems/ppzksnark/` covers both R1CS systems and GG aggregation.
- `test/systems/plonk/placeholder/` covers Placeholder with LPC/FRI and KZG.
- `test/commitment/` covers commitment lifecycle and failure cases.
- `test/transcript/`, `test/math/`, and `test/routing_algorithms/` cover the
  supporting modules.
