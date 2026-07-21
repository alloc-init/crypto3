# Poseidon1 and Sponge Construction {#hashes_poseidon1}

@tableofcontents

This page describes the standard Poseidon1 implementation and the shared sponge construction used by it.

## Implementation Map

The implementation is split across these headers:

|File|Purpose|
|----|-------|
|`detail/poseidon1/poseidon1_policy.hpp`|Poseidon1 policy, field/state/rate/digest definitions, and policy concept.|
|`detail/poseidon1/poseidon1_permutation.hpp`|Dense Poseidon1 permutation.|
|`detail/poseidon1/poseidon1_optimized_permutation.hpp`|Optimized Poseidon1 permutation equivalent to the dense one.|
|`detail/poseidon1/poseidon1_constants.hpp`|Dense constants adapter.|
|`detail/poseidon1/poseidon1_optimized_constants.hpp`|Optimized constants derived from dense constants.|
|`detail/poseidon_common/poseidon_sponge.hpp`|Shared sponge construction and sponge mode aliases.|

Tests for the policy, dense/optimized permutation equivalence, sponge modes, and digest output are in
`libs/hash/test/poseidon.cpp`.

## Poseidon1 Policy

The standard Poseidon1 policy is defined by `nil::crypto3::hashes::detail::poseidon1_policy`.

The policy parameters are:

|Parameter|Meaning|
|---------|-------|
|`FieldType`|Field used by the permutation and sponge. One sponge word is one field element.|
|`Security`|Security parameter used to select the checked-in Poseidon1 round count.|
|`Rate`|Number of field elements absorbed and squeezed per sponge block.|
|`Capacity`|Number of capacity field elements. Currently only capacity `1` has checked-in constants.|
|`DigestBits`|Requested digest size in bits. The policy rounds this up to whole field elements.|

The policy derives:

```cpp
word_bits = FieldType::value_bits
state_words = Rate + Capacity
block_words = Rate
digest_words = ceil(DigestBits / word_bits)
```

For compatibility, `DigestBits` defaults to `FieldType::value_bits`, so existing `poseidon1_policy<Field, Security,
Rate>` uses still return one field element.

If `digest_words == 1`, the digest type is `word_type`. If `digest_words > 1`, the digest type is
`std::array<word_type, digest_words>`.

## Output Limitation

The current sponge intentionally supports only:

```text
OUT <= RATE
```

Equivalently:

```cpp
digest_words <= block_words
```

`digest()` copies one squeezed rate block and converts its first `digest_words` elements into the digest. Supporting
larger outputs requires a multi-squeeze digest path.

Example over a 254-bit field:

```cpp
using policy = poseidon1_policy<field_type, 128, 2, 1, 256>;
```

This requests 256 digest bits. Since one field element has 254 bits, the policy derives:

```text
digest_words = 2
digest_type = std::array<word_type, 2>
```

## Poseidon1 Permutation

The standard dense permutation is `poseidon1_permutation<Policy>`.

It follows the original Poseidon1 round structure:

```text
RF / 2 full rounds
RP partial rounds
RF / 2 full rounds
```

Each dense round uses:

```text
AddRoundConstants -> S-box -> MDS
```

The optimized permutation is `poseidon1_optimized_permutation<Policy>`. It is algebraically equivalent to the dense
permutation, but rewrites the partial-round section so most partial rounds use sparse linear layers and scalar round
constants on state element `0`. The optimized constants are derived from the same dense MDS matrix and round constants.

## Shared Sponge Construction

The shared sponge is defined by `nil::crypto3::hashes::detail::poseidon_sponge_construction`.

It is parameterized by:

```cpp
template<
    typename PolicyType,
    typename PermutationType,
    poseidon_sponge_absorb_mode AbsorbMode = poseidon_sponge_absorb_mode::overwrite,
    poseidon_sponge_padding_mode PaddingMode = poseidon_sponge_padding_mode::pad10>
class poseidon_sponge_construction;
```

`PolicyType` supplies field/state/rate/digest types. `PermutationType` supplies:

```cpp
PermutationType::permute(state)
```

This keeps the sponge independent from the concrete permutation. The same construction can be used with dense
Poseidon1, optimized Poseidon1, and future Poseidon2 permutations.

## Absorb Modes

Two absorb modes are implemented:

|Mode|Rate update|
|----|-----------|
|`overwrite`|`state[i] = input[i]`|
|`add`|`state[i] += input[i]`|

The named production aliases currently use overwrite mode:

```cpp
poseidon_pad10_sponge_construction
poseidon_padding_free_sponge_construction
```

The generic template also supports add mode for callers that explicitly instantiate it.

## Padding Modes

Two padding modes are implemented.

### Pad10

`pad10` is the variable-length-safe mode.

If the final block is partial, the input words are absorbed into the rate, a sentinel `1` is placed at the next rate
cell, the remaining rate cells are zeroed in overwrite mode, and the state is permuted:

```text
[a] with RATE = 3 -> [a, 1, 0 | capacity...] -> P
```

If the final block is full, there is no unused rate cell for the sentinel. Instead, the first capacity cell is
incremented before the final permutation:

```text
[a, b] with RATE = 2 -> [a, b | capacity[0] + 1] -> P
```

The capacity cell is incremented rather than overwritten because it may already contain accumulated sponge state from
previous blocks.

### Padding-Free

`padding_free` does not add any length marker.

Full blocks are absorbed and permuted as-is. A partial final block absorbs only its present words and leaves the
remaining rate cells unchanged:

```text
state = P([a, b | capacity...])
final partial [c] -> [c, state[1] | capacity...] -> P
```

Padding-free mode is safe only when the input length is fixed by the caller's protocol context. It is not
collision-resistant for arbitrary variable-length input.

## Delayed Full Block

The streaming sponge delays the latest full block.

This is necessary because Pad10 treats an intermediate full block differently from a final full block:

```text
intermediate full block -> absorb rate, permute
final full block        -> absorb rate, increment first capacity cell, permute
```

`absorb(block)` therefore flushes the previous pending full block and stores the new block as pending. Finalization then
decides whether the pending block is final.

For input `[a, b, c, d]` with rate `2` and Pad10:

```text
absorb [a, b] as normal block -> P([a, b | cap])
final [c, d]                  -> P([c, d | cap + 1])
```

For padding-free mode, the same input is:

```text
absorb [a, b] as normal block -> P([a, b | cap])
final [c, d]                  -> P([c, d | cap])
```

## Squeezing and Digest

`squeeze()` finalizes the sponge if needed, copies one rate block from the state, permutes the state for possible
subsequent squeezes, and returns the copied block.

`digest()` calls `squeeze()` once and converts the first `digest_words` elements into the policy digest type.

Because `digest()` currently squeezes once, `digest_words <= block_words` is an intentional implementation limit.
