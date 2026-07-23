# Usage {#pubkey_usage_manual}

@tableofcontents

## BLS Signing

The range-based public-key algorithms accept byte containers directly:

```cpp
#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/pubkey/algorithm/sign.hpp>
#include <nil/crypto3/pubkey/algorithm/verify.hpp>
#include <nil/crypto3/pubkey/bls.hpp>

using namespace nil::crypto3;
using namespace nil::crypto3::algebra;
using namespace nil::crypto3::pubkey;

using curve_type = curves::bls12_381;
using scheme_type = bls<bls_default_public_params<>, bls_mss_ro_version,
                        bls_basic_scheme, curve_type>;
using private_key_type = private_key<scheme_type>;
using public_key_type = public_key<scheme_type>;
using scalar_type = typename private_key_type::private_key_type;
using signature_type = typename public_key_type::signature_type;

int main() {
    const std::string text = "hello world";
    const std::vector<std::uint8_t> message(text.begin(), text.end());

    private_key_type private_key(
        random_element<typename scalar_type::field_type>());
    const signature_type signature = sign(message, private_key);

    public_key_type &public_key = private_key;
    assert(verify(message, signature, public_key));
}
```

The maintained buildable version is `libs/pubkey/example/bls.cpp`.
