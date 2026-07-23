# Implementation {#block_ciphers_impl}

@tableofcontents

## Cipher Policies

Block cipher types such as `block::aes`, `block::rijndael`, and
`block::shacal2` define their key, block, and word types together with the
single-block encryption and decryption operations. Their public headers are
under `nil/crypto3/block`.

The policy types expose compile-time parameters such as block size, key size,
and round count. Constructing a cipher schedules its key; processing modes hold
that cipher while data is transformed.

## Algorithms

The public range and iterator algorithms are declared in
`block/algorithm/encrypt.hpp` and `block/algorithm/decrypt.hpp`. They construct
the appropriate encryption or decryption processing mode and feed input through
the block stream processor.

```cpp
#include <nil/crypto3/block/aes.hpp>
#include <nil/crypto3/block/algorithm/encrypt.hpp>

using namespace nil::crypto3;

auto ciphertext = encrypt<block::aes<128>>(input, key);
```

## Stream Processing

`detail::block_stream_processor` converts integral input values into the word
and block representation required by the selected cipher. Complete blocks are
forwarded to `block::accumulator_set`; the processing mode determines whether
the cipher encrypts or decrypts them.

## Result Values

Range-returning algorithm overloads produce a conversion wrapper declared in
`cipher_value.hpp`. It owns or references the accumulator state and can convert
the transformed data to a compatible output container. Output-iterator
overloads write directly to the caller's destination.
