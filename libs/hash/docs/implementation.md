# Implementation {#hashes_impl}

@tableofcontents

## Hash Policies

Hash policy types define the digest, word, block, compressor, and construction
used by an algorithm. Public policies are available under `nil/crypto3/hash`,
including SHA-2, SHA-3, BLAKE2, RIPEMD, Poseidon, and other families selected by
the corresponding CMake options.

Policies are stateless descriptions. Construction types implement padding and
finalization, while compressor types update the internal state for each message
block.

## Algorithms

The `nil::crypto3::hash` overloads in `hash/algorithm/hash.hpp` accept iterator
pairs, ranges, initializer lists, output iterators, or a caller-provided hash
accumulator.

```cpp
#include <string>

#include <nil/crypto3/hash/algorithm/hash.hpp>
#include <nil/crypto3/hash/sha2.hpp>

using namespace nil::crypto3;

std::string digest = hash<hashes::sha2<256>>(std::string("hello world"));
```

## Stream Processing

Each hash selects a stream processor that repacks integral input values into
the policy's word and block types. Complete blocks are submitted to
`hashes::accumulator_set`; finalization handles the remaining input according
to the selected construction.

## Result Values

Range-returning overloads produce the conversion wrapper declared in
`hash_value.hpp`. It can convert the digest to compatible containers while
preserving the hash policy's fixed digest size. Output-iterator overloads write
the digest directly to a destination range.
