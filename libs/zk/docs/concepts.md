# Concepts {#zk_concepts}

@tableofcontents

## R1CS

A Rank-1 Constraint System contains equations of the form

\f[
    \langle A, X \rangle \cdot \langle B, X \rangle = \langle C, X \rangle.
\f]

`snark::r1cs_constraint` stores the three linear combinations, and
`snark::r1cs_constraint_system` stores the constraints and assignment sizes.
Variable index `0` denotes the constant one and is not present in an assignment.
Primary variables begin at index `1`; auxiliary variables follow them. The
primary and auxiliary input vectors therefore contain exactly
`primary_input_size` and `auxiliary_input_size` elements respectively.

Use `is_valid()` to check variable indices and `is_satisfied()` to check an
assignment before generating a proof.

## PLONK and Placeholder

A `snark::plonk_constraint_system` describes gates, copy constraints, and
lookups. A `snark::plonk_assignment_table` holds witness, public-input,
constant, and selector columns, while `snark::plonk_table_description` records
their sizes and the usable-row range.

Placeholder divides preprocessing by data visibility:

- Public preprocessing handles the constraint system, public assignments,
  table description, and commitment parameters. Its `common_data` is required
  by both the prover and verifier.
- Private preprocessing handles witness assignments.
- The prover consumes both preprocessing results.
- The verifier consumes `common_data`, the proof, the circuit description, and
  a verifier-side commitment-scheme instance.

The prover and verifier must agree on the constraint system, table description,
commitment parameters, and transcript initialization.

## Commitment batches

Stateful commitment schemes use `commitments::polys_evaluator` to organize
polynomials and evaluation points by batch. The usual lifecycle is:

1. Add polynomials with `append_to_batch()` or `append_many_to_batch()`.
2. Commit to each batch with `commit()`.
3. Register evaluation points with `append_eval_point()`.
4. Produce an evaluation proof with `proof_eval()`.
5. Configure an equivalent verifier instance and call `verify_eval()`.

Batch identifiers and evaluation-point order are protocol inputs. They must be
identical on both sides. The maintained LPC tests demonstrate how the verifier
restores batch sizes from the proof before verification.

## Fiat-Shamir transcripts

`transcript::fiat_shamir_heuristic_sequential<HashType>` turns an interactive
protocol into a deterministic challenge sequence. Calling `operator()` absorbs
data; `challenge<FieldType>()` and `challenges<FieldType>()` derive challenges
and advance transcript state.

The prover and verifier must use the same initial state, absorb the same values
in the same order, and request challenges in the same order. Byte-oriented
hashes marshal field and curve elements before absorption. Field-oriented
hashes such as Poseidon require initialization and challenges over the policy's
field.

## Trusted setup

KZG and the pairing-based R1CS systems require structured reference strings.
Some tests and examples construct parameters from a known secret such as
`alpha`; that is convenient for deterministic testing but is not a secure
production setup ceremony. Applications are responsible for obtaining and
validating parameters generated with an appropriate process.
