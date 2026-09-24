# Modified SAP SNARK: protocol specification

## 1. Scope

This specification defines the standalone SNARK `modified_sap_snark`.
SAP means Square Arithmetic Program. Its three operations are `SNARK.Setup`,
`SNARK.Prove` and `SNARK.Verify`, defined by the relations, key structures and
equations below.

The protocol uses BN254, the modified squaring relation, one odd canonical
public scalar, evaluation order `A, C, H, Z`, and exact pairing normalization.

The construction uses a circuit-specific trusted setup and a Fiat-Shamir
challenge. It includes no zero-knowledge blinding and makes no zero-knowledge
claim. Unspecified protocol details are listed in section 9.

## 2. Relation, assignment and polynomial domain

All circuit coefficients, witness entries, polynomial coefficients and
evaluations are in `Fr`, BN254's scalar field of order `r`. The curve's base
field is `Fp`, of characteristic `p`. The two fields are distinct.

The native public input is one `Fr::value_type`, named `u`. Its canonical integer
representative must be odd: `0 <= canonical(u) < r` and `canonical(u) % 2 == 1`.
Parity refers to this integer, not an arbitrary representative modulo `r`.

Assignment contract:

- `n >= 1` is the total number of witness entries, including slot zero.
- `witness` contains exactly `n` scalar values, with `witness[0] == u`.
- A linear-combination term with index `i` reads `witness[i]`, for `0 <= i < n`.
  There is no implicit constant-one entry.
- Use Crypto3's `assignment_layout::explicit_constant` for this direct indexing.
  Its name does not impose the value one. A term at index zero, including one
  made by a scalar-valued linear-combination constructor, multiplies `u`.

The constraint system stores ordered sparse row pairs `(a_j, c_j)` with:

```text
A_j(w) = sum_i a[j,i] * w[i]
C_j(w) = sum_i c[j,i] * w[i]
A_j(w)^2 - C_j(w) = u
```

These rows define the modified squaring constraint system. The first logical
row must be the public-input binding row: `A_0 = 0`, `C_0 = -w[0]`.
It enforces `0^2 - (-w[0]) = u`. Setup validates this row; it does not silently
insert variables or reinterpret a conventional R1CS assignment. A host-side
`witness[0] == u` check alone would not put that binding into the proved relation.

Let `q >= 1` count logical rows, including the binding row. Choose:

```text
m = smallest supported power of two >= max(2, q)
omega = math::unity_root<Fr>(m)
D = (1, omega, ..., omega^(m-1))
Z(X) = product_j (X - D[j]) = X^m - 1
```

Use `math::basic_radix2_domain<Fr>` and its natural domain-element order. Reject
unsupported sizes or size arithmetic overflow. Pad from `q` to `m` rows with
copies of the binding row; an all-zero row would incorrectly require `u == 0`.

Interpolate the padded coefficient columns so that `A^i(D[j]) = a[j,i]` and
`C^i(D[j]) = c[j,i]`. The witness determines:

```text
A(X) = sum_i w[i] * A^i(X)
C(X) = sum_i w[i] * C^i(X)
H(X) = (A(X)^2 - C(X) - u) / Z(X)
A(X)^2 - C(X) = Z(X)*H(X) + u
```

The numerator is divisible by `Z` exactly when all padded rows are satisfied.
Proving computes `H` internally and rejects a nonzero remainder. Coefficient
vectors use ascending powers: element `i` is the coefficient of `X^i`.

| Polynomial | Degree bound | Maximum coefficient count | Opening quotient coefficient count |
| --- | --- | --- | --- |
| `A`, `C` | `m - 1` | `m` each | `m - 1` each |
| `H` | `m - 2` | `m - 1` | `m - 2` |
| `Z` | exactly `m` | `m + 1` | `m` |

Here an opening quotient is `(U(X) - U(z)) / (X - z)`, (`U= A, C, H, Z`). Zero polynomials are
allowed. At `m = 2`, `H` is constant and its opening quotient is zero, requiring
no opening bases. Bounds describe logical sizes, independently of polynomial
containers trimming trailing zeros.

Conversion from ordinary R1CS must explicitly remap the original constant and
public-input indices and introduce any required auxiliaries. Recover
the constant-one value using a constrained binary decomposition:

```text
b[i] in {0, 1}, for 0 <= i < bit_length(r)
sum_i 2^i * b[i] == canonical(u), as integers with the sum < r
constant_one = b[0]
```

The verifier's oddness check then guarantees `b[0] == 1`. The circuit must enforce
both the field reconstruction and the integer range; field equality alone would
permit a noncanonical representative with different parity.

### 2.1 Ordinary R1CS frontend

The frontend is `snark::reductions::r1cs_to_modified_sap<FieldType>`, under
`snark/reductions/r1cs_to_modified_sap.hpp`. For this protocol, `FieldType` is
`Fr`. It accepts Crypto3's existing `r1cs_constraint_system<FieldType>` with
exactly one public input. Its interface contract is:

```cpp
using reduction = snark::reductions::r1cs_to_modified_sap<Fr>;
auto constraint_system = reduction::instance_map(r1cs);
auto witness = reduction::witness_map(r1cs, primary_input, auxiliary_input);
```

`instance_map` returns an owned `modified_sap_constraint_system<FieldType>`.
`witness_map` takes the existing `r1cs_primary_input<FieldType>` and
`r1cs_auxiliary_input<FieldType>` vectors and returns
`std::vector<typename FieldType::value_type>`, containing the complete modified
SAP witness. Both operations borrow their inputs as const references, preserve
them, and retain no references to them.

The source public-input vector contains exactly one element, `u`; the native
SNARK operations receive that element as a scalar. Circuit conversion depends
only on the source constraints, declared dimensions, and field modulus. The same
converted circuit and setup serve every satisfying assignment with an admissible
odd `u`. Polynomial construction and domain padding use the reduction defined
earlier in this section.

#### Source indexing and witness layout

Let `N = 1 + r1cs.auxiliary_input_size` be the number of source variables,
excluding its implicit constant, and let `M` be its number of constraints.
Source index zero means one, index one is `u`, and indices `2..N` are the source
auxiliaries in their original order. Let:

```text
k = bit_length(r)
h = number of 1 bits of r, excluding its most and least significant bits
g = number of maximal consecutive runs of 0 bits of r
s = h + g + 1
```

Since `r` is odd, its most and least significant bits are both one. The constant
recovery below uses `k` bits, `h` new prefix markers and `s` multiplication
constraints. Each multiplication introduces one additional square-conversion
variable. All witness ranges below are half-open:

| Witness range | Contents and order |
| --- | --- |
| `[0, N)` | Source assignment: `u`, then the original auxiliaries |
| `[N, N + k)` | Bits `b[0]..b[k-1]` of `canonical(u)`, least significant first |
| `[N + k, N + k + h)` | New comparison prefix markers, in most-to-least-significant traversal order |
| `[N + k + h, N + k + h + s)` | Square-conversion variables for the comparison products, in emission order |
| `[N + k + h + s, N + k + h + s + M)` | Square-conversion variables for the source constraints, in source row order |

Remap each source term without changing its coefficient:

```text
source index 0       -> witness index N       // recovered one: b[0]
source index i >= 1  -> witness index i - 1
```

In particular, a source constant coefficient becomes a coefficient of `w[N]`,
not `w[0]`. Scalar coefficients in linear combinations remain fixed field
elements; they do not require a witness variable of their own.

#### Multiplication and square conversion

Reuse the R1CS-to-SAP conversion for each multiplication `L * R = O`, where
`L`, `R` and `O` are the left, right and output linear combinations over the
explicit witness indices, respectively.
Introduce one variable `t`, whose honest value is `(L(w) - R(w))^2`, and emit
these two modified SAP rows, in this order:

```text
a = L + R,    c = 4*O + t - w[0]
a = L - R,    c = t - w[0]
```

With the binding row enforcing `w[0] = u`, these rows state
`(L + R)^2 = 4*O + t` and `(L - R)^2 = t`. Subtraction gives `4*L*R = 4*O`,
which is equivalent to the source multiplication because `r` is odd. Use this
same conversion for comparison products and source constraints. It adds two
rows and one variable for each product, including products with empty or
constant linear combinations.

More generally, an ordinary square equation `a(w)^2 = d(w)` becomes the modified
row `(a, d - w[0])`. The frontend returns canonical sparse rows: combine
duplicate indices, remove zero coefficients, and sort terms, preserving every
row and its position. Only the modified SAP binding and constant-recovery rows
are added to the converted source rows; ordinary SAP's additional public-input
independence rows belong to that separate reduction.

#### Canonical constant recovery

Write `r[j]` for bit `j` of the integer modulus. First constrain each witness bit
to be Boolean by emitting `(a, c) = (b[j], b[j] - w[0])` for `j = 0..k-1`.
Then emit one reconstruction row:

```text
B = sum_(j=0..k-1) 2^j * b[j]
a = 0,    c = -B
```

This row enforces `B = u` in Fr. Its coefficients `2^j` are field elements.
The integer represented by these bits must additionally satisfy `B < r`.
Enforce this bound by comparing from the most significant bit downward:

1. Start with marker `e = b[k-1]`, reusing that bit's witness slot. For a prefix
   that has passed the comparison, `e = 1` means it equals the corresponding
   prefix of `r`; `e = 0` means it is already smaller.
2. Traverse positions `k-2` down to `1`. At each maximal run of zero bits of
   `r`, emit `e * sum_(j in run) b[j] = 0`, leaving `e` unchanged. At each one
   bit, allocate a new marker `e_next`, emit `e * b[j] = e_next`, and replace
   `e` with `e_next`.
3. At position zero, emit `e * b[0] = 0`. Because `r[0] = 1`, this requires a
   strict decrease at that position or an earlier one, excluding equality.

Emit the two square-conversion rows for each product immediately in traversal
order. Every marker is Boolean by induction from the input bits and the product
equations. Each zero run contains fewer than `r` bits, so a sum of Boolean ones
cannot wrap to zero in Fr. Thus a zero-run check with `e = 1` forces every bit in
that run to zero; with `e = 0`, lower bits are unrestricted by the comparison.

These constraints require no constant-one witness: they use only bit variables,
prefix markers, zero linear combinations, fixed coefficients, and the binding
to `u`. Booleanity, reconstruction and the bound together force the bits to be
the canonical representation of `u`. Its required oddness then establishes
`b[0] = 1`, making the source constant remapping sound. In particular, bits of
`u + r` cannot substitute for those of `u`, even when both fit in `k` bits and
give the same field value.

The witness map computes the bits from `u.to_integral()`, each new prefix marker
as the product of the preceding marker and the corresponding bit, and each
square-conversion variable as `(L(w) - R(w))^2`. It then computes the source
square-conversion variables in source row order. These computations provide
witness values; the emitted equations independently enforce their correctness.

#### Row order and dimensions

The logical rows are ordered as follows:

| Row range | Contents |
| --- | --- |
| `[0, 1)` | Public-input binding: `(a, c) = (0, -w[0])` |
| `[1, k + 1)` | Booleanity of `b[0]..b[k-1]` |
| `[k + 1, k + 2)` | Reconstruction of `u` from the bits |
| `[k + 2, k + 2 + 2*s)` | Comparison products, two rows per product |
| `[k + 2 + 2*s, k + 2 + 2*s + 2*M)` | Source products, two rows per original constraint |

Consequently:

```text
n = N + k + h + s + M
q = k + 2 + 2*s + 2*M
```

For BN254 Fr, `k = 254`, `h = 99`, `g = 53` and `s = 153`. The frontend
therefore produces `n = N + 506 + M` witness entries and `q = 562 + 2*M` logical
rows. These counts exclude domain padding. The radix-two domain size `m` is
determined from `q` as above. Even a source with no multiplication constraints
retains the binding and constant-recovery rows.

#### Validation and failure behavior

Both maps throw `std::invalid_argument` for a source public-input count other
than one, out-of-range source indices, size arithmetic or container-size
overflow, or dimensions unsupported by the modified SAP domain. Check every
raw source index against `0..N`, including zero-coefficient and cancelling
terms, before normalization or evaluation. Unsorted and repeated terms with
valid indices are accepted and normalized using the existing linear-combination
operations. Empty combinations and unused source variables are allowed.

The witness map additionally requires `primary_input.size() == 1`,
`auxiliary_input.size() == r1cs.auxiliary_input_size`, an odd canonical `u`, and
satisfaction of every source row under the ordinary implicit-one convention.
Otherwise it throws `std::invalid_argument` and returns no witness. Check all
dimensions and indices before invoking source evaluation. Validation applies in
release builds as well as debug builds. Allocation failures propagate.

For a structurally valid source and admissible `u`, every satisfying source
assignment has the extension described above. Conversely, any assignment
satisfying the converted relation has `w[0] = u` and `w[N] = 1`; projecting its
first `N` entries recovers a satisfying source assignment. This converse applies
to arbitrary supplied witnesses, independently of the witness map.

#### Worked example

Take the source constraint `(x + 1) * y = u`, with public input `[3]` and
auxiliaries `[2, 1]`. Source indices are `0 = one`, `1 = u`, `2 = x`, `3 = y`,
so `N = 3` and `M = 1`. The modified witness begins with `[3, 2, 1]`, followed
by the bits of 3: `b[0] = w[3] = 1`, `b[1] = w[4] = 1`, and all other bits zero.

The source row becomes `L = w[1] + w[3]`, `R = w[2]`, `O = w[0]`.
Its square-conversion variable is `t = (3 - 1)^2 = 4`. The two rows evaluate to:

```text
(3 + 1)^2 - (4*3 + 4 - 3) = 16 - 13 = 3
(3 - 1)^2 - (4 - 3)       =  4 -  1 = 3
```

Over BN254 Fr, `t` occupies witness slot 509, the total witness size is 510,
the logical row count is 564, and the padded domain size is 1024. Changing the
public input to 5 and the auxiliaries to `[4, 1]` uses the same converted circuit;
only the witness values change.

## 3. Groups and pairing convention

Use `algebra::curves::alt_bn128_254` and its `g1_type<>`, `g2_type<>`, scalar,
base and target-field types. Fix Crypto3's standard G1 and G2 generators.

```text
[x]_1 = x * G1_generator
[x]_2 = x * G2_generator
gT = e_exact(G1_generator, G2_generator)
[x]_T = gT^x
E = (p^12 - 1) / r
e_exact(P, R) = MillerLoop(P, R)^E
```

G1 and G2 use additive notation. GT uses multiplication: `[0]_T` is the field
element one, while the field element zero is not a GT member. All GT equations
below use multiplication and exponentiation.

Normalize Crypto3's optimized BN254 final exponent using:

```text
t = 0x44E992B44A6909F1
c = 2*t*(6*t^2 + 3*t + 1)
exact_result = optimized_result ^ inverse(c mod r, r)
```

Select this through an explicit compile-time pairing policy. Setup
and verification must use the same convention, including `gT` and every GT key
element.

For nonzero Miller-loop values, normalization must equal direct exponentiation
by `E`. The final-exponentiation operation has an optional-result
contract: zero input fails, and input one returns one. Pairing with a valid
identity point must give the GT identity; identity handling cannot depend on an
affine conversion that divides by zero.

## 4. Setup and key contents

`SNARK.Setup` takes the circuit, not a particular witness or value of `u`.
It samples independent secret scalars `tau`, `gamma`, and
`alpha_A, alpha_C, alpha_H, alpha_Z` from Fr. Resample `tau` until it is nonzero
and `Z(tau) != 0`, and resample `gamma` until it is nonzero. Additional rejection
rules for degenerate samples are unspecified; see section 9.

These checks prevent a degenerate reference string. If `tau == 0`, every positive
power of `tau` vanishes, most polynomial-query entries collapse to the group
identity, and `tau_g2` exposes the degeneration. If `Z(tau) == 0`, `tau` lies in
the SAP evaluation domain and `alpha_z_vanishing_gt` becomes the GT identity.
A zero `gamma` has no inverse, so `gamma_inverse_g2` could not be constructed.

The supplied randomness source is consumed by reference and is not retained.
Production callers must supply cryptographically secure randomness; deterministic
sources are for reproducible checks. The source must satisfy the uniform random
bit generator requirements. Setup samples canonical integers uniformly from
`[0, r - 1]` in `tau, gamma, alpha_A, alpha_C, alpha_H, alpha_Z` order. It finishes
resampling `tau` before sampling `gamma`, then samples the four alphas. Keys contain
group encodings, never the secret scalars themselves.

Setup constructs the following G1 query vectors for the proving key:

| Field | Entry | Length / index range |
| --- | --- | --- |
| `W` | `[gamma*(alpha_A*A^i(tau) + alpha_C*C^i(tau))]_1` | `n`, `0 <= i < n` |
| `H_query` | `[gamma*alpha_H*tau^i]_1` | `m - 1`, `0 <= i < m - 1` |
| `T_A`, `T_C` | `[alpha_U*tau^i]_1` | `m - 1` each, `0 <= i < m - 1` |
| `T_H` | `[alpha_H*tau^i]_1` | `m - 2`, `0 <= i < m - 2` |
| `T_Z` | `[alpha_Z*tau^i]_1` | `m`, `0 <= i < m` |

The verification key contains three G2 values and five GT values:

| Field | Value |
| --- | --- |
| `g2_one` | `[1]_2`, the fixed generator |
| `tau_g2` | `[tau]_2` |
| `gamma_inverse_g2` | `[gamma^-1]_2` |
| `alpha_z_vanishing_gt` | `[alpha_Z*Z(tau)]_T` |
| `alpha_gt`, ordered `A, C, H, Z` | `[alpha_U]_T` for each `U` |

Public metadata is `num_variables = n`, `domain_size = m`, and
`circuit_digest`, one `Fp` value for the transcript profile. These are
binding data stored alongside the group elements defined above.
The digest binds the circuit, including logical row count, coefficients, witness
indexing and domain convention; it contains no particular witness or public input.

The proving key owns a canonical copy of the logical constraint system, the query
vectors, and a copy of the verification key. It can reconstruct padding and the
domain from this data. Canonical sparse rows combine duplicate indices, omit zero
coefficients and sort indices without changing row order. The canonical binding
row is checked after that normalization.

The keypair is `std::pair<proving_key_type, verification_key_type>`. Its second
member equals the verification-key copy in the first member. Returning from setup
must not leave references to the caller's circuit or to temporary storage.
Calls to prove and verify borrow their inputs as const references and do not
retain witness data. FFT and pairing caches are derived data, not additional
protocol parameters.

## 5. Proof and verification equations

The proof owns `P` and `Q` in G1, followed by the four Fr evaluations
`v_A, v_C, v_H, v_Z`. The public input `u` is a separate argument: the verifier
receives the expected value from the caller and uses it in every check and in
the transcript.

After deriving `A, C, H, Z`, `SNARK.Prove` computes:

```text
P = sum_(i=0..n-1) w[i]*W[i]
    + sum_(i=0..m-2) h[i]*H_query[i]
z = Challenge(public_parameters, P, u)
v_U = U(z), for U in (A, C, H, Z)
Q_U(X) = (U(X) - v_U) / (X - z)
Q = sum_U sum_i coefficient(Q_U, i)*T_U[i]
```

Pad coefficient vectors with zero only as needed for their corresponding query
lengths. Empty sums give the G1 identity. Quotient construction is polynomial
division, so `z == 0` and `z` belonging to the circuit domain are valid; it must
not divide by `z` or `Z(z)`.

After validation, `SNARK.Verify` recomputes the same challenge and requires:

```text
canonical(u) % 2 == 1
v_A^2 - v_C == v_H*v_Z + u

e_exact(P, gamma_inverse_g2) * alpha_z_vanishing_gt
  == e_exact(Q, tau_g2)
     * e_exact(Q, g2_one)^(-z)
     * product_(U in A,C,H,Z) alpha_gt[U]^v_U
```

The negative exponent is in the order-`r` target group. Opening `Z` is required
by this verification equation even though native code could evaluate `X^m - 1`.

Validate curve and order-`r` subgroup membership of proof points and verification
key points before pairing. GT constants must be nonzero and satisfy `x^r == 1`.
The stored G2 generator must equal the fixed generator, and `tau_g2` and
`gamma_inverse_g2` must be nonidentity. Do not blanket-reject identity proof
points: valid commitments can cancel. Group-order checks must use integer
multiplication by `r`, not a conversion to an Fr scalar that would reduce `r`
to zero.

Structural and membership checks do not certify an honestly generated CRS.
The caller selects the trusted circuit-specific verification key; verifying a
proof under a key supplied by the prover alone does not establish that binding.

## 6. Native API

The scheme facade is `snark::modified_sap_snark<Policy>`, under
`snark/systems/ppsnark/modified_sap_snark.hpp`. The compile-time policy identifies
the curve, exact pairing convention and transcript profile. There is no runtime
switch between transcript conventions. The interface contract is:

```cpp
using scheme = snark::modified_sap_snark<policy>;
auto keys = scheme::generate(constraint_system, random_source);
auto proof = scheme::prove(keys.first, u, witness);
bool valid = scheme::verify(keys.second, u, proof);
```

The scheme exposes these aliases:

| Alias | Contract |
| --- | --- |
| `constraint_system_type` | The modified-SAP constraint system over Fr |
| `primary_input_type` | One Fr value, `u` |
| `auxiliary_input_type` | A vector of Fr values: the full explicitly indexed witness, including `w[0]` |
| `proof_type` | The two G1 points and four Fr evaluations in section 5 |
| `proving_key_type`, `verification_key_type` | The owned data in section 4 |
| `keypair_type` | `std::pair<proving_key_type, verification_key_type>` |

The full-witness meaning of `auxiliary_input_type` is specific to this scheme.
The generic adapter contract is `zk::generate<scheme>(cs, rng)`,
`zk::prove<scheme>(pk, u, witness)` and `zk::verify<scheme>(vk, u, proof)`, each
forwarding to the corresponding scheme operation. Setup requires an explicit
randomness source.

Proving is deterministic for a fixed key, public input and witness. It consumes
no additional randomness and leaves those inputs unchanged. A complete proof is
returned only after both commitment and opening computation succeed.

Failure behavior:

| Operation | Failure behavior |
| --- | --- |
| `generate` | Throw `std::invalid_argument` for invalid indices, missing binding row, unsupported dimensions or inconsistent circuit structure. Reject size overflow before allocation. |
| `prove` | Throw `std::invalid_argument` for inconsistent key dimensions, wrong witness length, even `u`, `w[0] != u`, or an unsatisfied relation. Return no partial proof. |
| `verify` | Return `false` for invalid public input, malformed key/proof, failed membership checks, pairing failure or either failed verification equation. |

Allocation and randomness-source failures propagate; they are not converted
into proof rejection. Validation is not conditional on assertions or debug mode.
Native field values already denote residues; decoding must reject an
out-of-range integer before constructing a field value and losing that evidence.

### 6.1 Ordinary R1CS example

Supply a Crypto3 ChaCha generator seeded from cryptographic entropy. This example
converts `(x + 1) * y = u` and runs setup once, then proves the assignments
`(u, x, y) = (3, 2, 1)` and `(5, 4, 1)` with the same keys. The frontend accepts
the ordinary public and private input vectors; the SNARK receives the scalar
`u` and the complete converted witness.

```cpp
#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/math/linear_variable.hpp>
#include <nil/crypto3/random/chacha_urbg.hpp>
#include <nil/crypto3/zk/snark/reductions/r1cs_to_modified_sap.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/transcript.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap_snark.hpp>

bool prove_r1cs_assignments(nil::crypto3::random::chacha_urbg<> &random_source) {
    namespace snark = nil::crypto3::zk::snark;
    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using field_type = curve_type::scalar_field_type;
    using value_type = field_type::value_type;
    using variable_type = nil::crypto3::math::linear_variable<field_type>;
    using policy = snark::modified_sap_policy<curve_type, snark::modified_sap_bn254_exact_pairing_policy,
                                              snark::modified_sap_bn254_poseidon_transcript_policy>;
    using scheme = snark::modified_sap_snark<policy>;
    using frontend = snark::reductions::r1cs_to_modified_sap<field_type>;

    snark::r1cs_constraint_system<field_type> source;
    source.primary_input_size = 1;
    source.auxiliary_input_size = 2;
    // Source indices: 0 = one, 1 = u, 2 = x, 3 = y.
    snark::r1cs_constraint<field_type> row {variable_type(0), variable_type(3), variable_type(1)};
    row.a.add_term(variable_type(2));    // (1 + x)*y = u.
    source.add_constraint(row);

    const auto constraint_system = frontend::instance_map(source);
    const auto keys = scheme::generate(constraint_system, random_source);
    for (const auto &u : {value_type(3), value_type(5)}) {
        const auto witness = frontend::witness_map(source, {u}, {u - value_type::one(), value_type::one()});
        const auto proof = scheme::prove(keys.first, u, witness);
        if (!scheme::verify(keys.second, u, proof)) {
            return false;
        }
    }
    return true;
}
```

## 7. Transcript requirements

### 7.1 Poseidon profile

Both the circuit digest and proof challenge use Crypto3 Poseidon1 over BN254's
base field Fp. The exact Crypto3 type is:

```cpp
hashes::poseidon1<
    hashes::detail::poseidon1_policy<
        algebra::fields::alt_bn128_base_field<254>, 128, 2>>
```

This selects a width-three state, rate two, capacity one, an `x^5` S-box, eight
full rounds and 56 partial rounds. The round schedule is four full rounds, 56
partial rounds, then four full rounds. Each ordinary round adds its round
constants, applies the S-box to every cell in a full round or only cell zero in
a partial round, then multiplies by the checked-in dense MDS matrix. The selected
`hashes::poseidon1` type executes Crypto3's optimized partial-round permutation.
It is algebraically equivalent to the dense schedule described above and produces
the same digest using the same checked-in BN254 base-field constants. The dense
schedule defines the canonical result; implementations may use either equivalent
permutation.

Hash one complete field-element sequence as one message. Initialize the state to
`[0, 0, 0]`, use overwrite absorption, and apply Crypto3's `pad10` mode. Each
non-final pair replaces state cells zero and one and is followed by a permutation.
For a final one-element block, overwrite cell zero with the element and cell one
with `1`, then permute. For a final full pair, overwrite cells zero and one, add
`1` to capacity cell two, then permute. An empty message sets cells zero and one
to `1` and `0` before the permutation. The digest is cell zero after the final
permutation. No transcript operation performs a separate hash per absorbed item.

### 7.2 Tags and primitive encodings

`tag(s)` is the nonnegative integer represented by the ASCII bytes of `s` in
big-endian order, embedded into Fp. Every tag below is shorter than the Fp
modulus width.

| Meaning | String | Canonical integer |
| --- | --- | --- |
| Protocol | `modified-sap-snark` | `0x6d6f6469666965642d7361702d736e61726b` |
| Version | `v1` | `0x7631` |
| Curve | `bn254` | `0x626e323534` |
| Poseidon profile | `poseidon1-fp-128-r2-c1` | `0x706f736569646f6e312d66702d3132382d72322d6331` |
| Pairing convention | `pairing-exact-e` | `0x70616972696e672d65786163742d65` |
| Circuit digest | `circuit-digest` | `0x636972637569742d646967657374` |
| Proof challenge | `proof-challenge` | `0x70726f6f662d6368616c6c656e6765` |
| Circuit row | `row` | `0x726f77` |
| A combination | `A` | `0x41` |
| C combination | `C` | `0x43` |

Encode a nonnegative size or index as one Fp element with the same integer value.
Reject values greater than or equal to `p`; native `size_t` values are always
smaller for BN254. Encode an Fr element as its canonical integer in `[0, r - 1]`,
embedded unchanged into Fp. In particular, the Fr value `-1` is encoded as
`r - 1`, rather than `p - 1`. BN254 has `r < p`, so this embedding is injective.

Encode a G1 point as three Fp elements. The identity is `[0, 0, 0]`. For a
nonidentity point, convert it to affine coordinates and encode `[1, x, y]`.

Encode a G2 point as five Fp elements. The identity is `[0, 0, 0, 0, 0]`. For a
nonidentity point with affine Fp2 coordinates, encode:

```text
[1, x.c0, x.c1, y.c0, y.c1]
```

Here `c0` and `c1` are the extension components in Crypto3 `data[0]`, `data[1]`
order. Test identity before affine conversion so identity encoding never performs
an inversion.

Encode a GT element as its 12 Fp components in Crypto3 tower order:

```text
for outer in 0..1:
    for middle in 0..2:
        for inner in 0..1:
            value.data[outer].data[middle].data[inner]
```

This is a fixed-width encoding and therefore needs no identity flag. Verification
rejects a zero GT value before deriving the challenge.

### 7.3 Circuit digest

Normalize the logical constraint system before hashing. Preserve logical row
order; within each A or C combination, combine duplicate indices, remove zero
coefficients and sort terms by increasing witness index. Let `n` be the witness
size including slot zero, `logical_rows` the logical row count, and `m` the padded
radix-two domain size. Hash the following single Fp sequence:

```text
tag("modified-sap-snark")
tag("v1")
tag("bn254")
tag("poseidon1-fp-128-r2-c1")
tag("circuit-digest")
n
logical_rows
m

for each logical row j in source order:
    tag("row")
    j
    tag("A")
    number of A terms
    for each A term:
        witness index
        coefficient encoded from Fr into Fp
    tag("C")
    number of C terms
    for each C term:
        witness index
        coefficient encoded from Fr into Fp
```

The resulting Fp digest is stored in the verification key. Padding rows are not
absorbed individually: their convention is fixed by this protocol version, and
`logical_rows` and `m` bind their number. No witness or public input is part of
this digest.

### 7.4 Proof challenge

Validate the proof and verification-key group elements before encoding them.
Hash the following single Fp sequence in exactly this order:

```text
tag("modified-sap-snark")
tag("v1")
tag("bn254")
tag("poseidon1-fp-128-r2-c1")
tag("pairing-exact-e")
tag("proof-challenge")
circuit_digest
num_variables
domain_size
3
encode_G2(g2_one)
encode_G2(tau_g2)
encode_G2(gamma_inverse_g2)
5
encode_GT(alpha_z_vanishing_gt)
encode_GT(alpha_gt[A])
encode_GT(alpha_gt[C])
encode_GT(alpha_gt[H])
encode_GT(alpha_gt[Z])
encode_G1(P)
u encoded from Fr into Fp
```

The literal `3` and `5` frame the G2 and GT lists. The challenge does not absorb
`Q` or the claimed evaluations because they are computed after `z`. Every actual
verification-key field is absorbed directly; there is no cached
verification-parameter digest.

Let `h` be the canonical integer representative of the resulting Fp digest.
The verifier and prover derive the challenge as `z = Fr(h mod r)`. This uses the
full Fp output. It is reduction rather than bit truncation or rejection sampling.

## 8. Correctness

For a satisfying assignment, the polynomial identity in section 2 holds for all
`X`. Evaluating it at the transcript challenge `z` gives the verifier's arithmetic
equation. The binding row and the canonical odd public input provide its public
input conditions.

For each `U` in `(A, C, H, Z)`, polynomial division gives:

```text
U(tau) = (tau - z)*Q_U(tau) + v_U
```

By bilinearity and the key definitions, the left side of the pairing equation
is `gT` raised to `sum_U alpha_U*U(tau)`. The right side is `gT` raised to:

```text
(tau - z)*sum_U alpha_U*Q_U(tau) + sum_U alpha_U*v_U
```

The exponents are equal by the quotient identities, so an honestly generated
proof satisfies the pairing equation. This argument also covers `tau == z`
without dividing by `tau - z`. It establishes algebraic correctness, not a
knowledge-soundness proof.

## 9. Unspecified protocol details

| Area | Unspecified details |
| --- | --- |
| Setup sampling | Rejection rules beyond `tau != 0`, `Z(tau) != 0`, and `gamma != 0`. |
| Serialization | Wire format, version identifiers, byte order, canonical decoding and resource limits. |
