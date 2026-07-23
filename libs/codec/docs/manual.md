# Manual # {#codec_manual}

@tableofcontents

## Basic Example

```cpp
#include <string>

#include <nil/crypto3/codec/base.hpp>
#include <nil/crypto3/codec/algorithm/encode.hpp>

using namespace nil::crypto3::codec;

int main(int argc, char *argv[]) {
    std::string data = "hello world";
    
    std::string result = encode<base64>(data);
    
    return !(result == "aGVsbG8gd29ybGQ=");
}
```
