# Introduction {#zk_introduction}

Crypto3.ZK is a header-only C++ library for constructing and verifying
zero-knowledge proofs. Public headers are under `nil/crypto3/zk`, primarily in
the `nil::crypto3::zk` and `nil::crypto3::zk::snark` namespaces.

@tableofcontents

## Proof systems

The library contains the following proof-system implementations:

- `snark::r1cs_ppzksnark`, a pairing-based preprocessing zkSNARK for R1CS.
- `snark::r1cs_gg_ppzksnark`, a Groth16-style R1CS preprocessing zkSNARK. It
  supports basic proving, proof aggregation, and encrypted-input modes through
  `snark::proving_mode`.
- PLONK Placeholder, exposed as public/private preprocessors, a prover, and a
  verifier. Placeholder includes permutation, gate, lookup, quotient, and dFRI
  argument components and supports LPC/FRI and KZG commitment backends.

## Arithmetization

The R1CS data model consists of linear combinations, constraints, a constraint
system, primary input, and auxiliary input. R1CS instances can be reduced to:

- QAP, a Quadratic Arithmetic Program.
- SAP, a Square Arithmetic Program.

The PLONK data model includes variables, constraints, gates, copy constraints,
lookup gates and tables, assignment tables, and polynomial tables.

## Polynomial commitments

Implemented polynomial commitment primitives include:

- FRI and the FRI-backed List Polynomial Commitment (LPC).
- KZG, batched KZG, KZG v2, and KZG IPP2 support.
- Pedersen commitments.
- Powers of Tau, proof-of-knowledge, proof-of-work, knowledge-commitment, and
  R1CS GG-ppzkSNARK MPC support.

## Supporting modules

Crypto3.ZK also provides sequential Fiat-Shamir transcripts, expression trees
and cached DAG evaluation, generic `generate`, `prove`, `verify`, and
`aggregate` adapters, and Benes and AS-Waksman routing networks.
