# Usage {#hashes_usage_manual}

@tableofcontents

## Quick Start

The range overload of `nil::crypto3::hash` hashes a container directly and
returns a value convertible to a compatible output container:

```cpp
#include <iostream>
#include <string>

#include <nil/crypto3/hash/algorithm/hash.hpp>
#include <nil/crypto3/hash/sha2.hpp>

int main() {
    const std::string input = "A quick brown fox jumps over the lazy dog";
    const std::string digest =
        nil::crypto3::hash<nil::crypto3::hashes::sha2<256>>(input);

    std::cout << digest << '\n';
}
```

Iterator overloads are available when only part of a range should be hashed:

```cpp
const std::string digest =
    nil::crypto3::hash<nil::crypto3::hashes::sha2<256>>(
        input.begin(), input.end());
```

See `libs/hash/example/hashes.cpp` for maintained SHA-2, SHA-3, Keccak, and
Poseidon examples.
