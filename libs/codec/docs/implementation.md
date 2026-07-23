# Implementation {#codec_impl}

@tableofcontents

## Codecs

Codec types describe both encoding and decoding modes. The maintained codec
families are declared in:

- `nil/crypto3/codec/base.hpp` for Base32, Base58, and Base64
- `nil/crypto3/codec/hex.hpp` for hexadecimal encoding

Each codec exposes `stream_encoder_type` and `stream_decoder_type`. These modes
bind the codec policy to the stream processor used by the generic algorithms.

## Algorithms

The public entry points are `nil::crypto3::encode` and
`nil::crypto3::decode`, declared in `codec/algorithm/encode.hpp` and
`codec/algorithm/decode.hpp`. Both support iterator pairs, ranges, initializer
lists, output iterators, and caller-provided codec accumulators.

```cpp
#include <string>

#include <nil/crypto3/codec/algorithm/decode.hpp>
#include <nil/crypto3/codec/algorithm/encode.hpp>
#include <nil/crypto3/codec/base.hpp>

using namespace nil::crypto3;
using namespace nil::crypto3::codec;

std::string encoded = encode<base64>(std::string("hello world"));
std::string decoded = decode<base64>(encoded);
```

Range-returning overloads produce a conversion wrapper from `codec_value.hpp`.
The wrapper can be converted to a compatible output container. Output-iterator
overloads write directly to a destination range.

## Streaming

Fixed-size codecs use `fixed_block_stream_processor`; codecs whose output size
depends on the full input use `varlength_block_stream_processor`. Processors
feed `codec::accumulator_set`, which stores the transformed blocks and handles
codec finalization.

The stateful interfaces in `codec_state.hpp` expose the same accumulator-based
machinery for callers that need to supply input incrementally.
