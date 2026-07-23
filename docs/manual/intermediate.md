\page BLS BLS Signing

# BLS Signing

This guide follows the maintained example in `libs/pubkey/example/bls.cpp`. It
creates a BLS12-381 private key, signs a byte range, and verifies the signature.

Weighted threshold BLS modes described by older versions of this page are not
part of the current public-key API. Weighted Shamir secret sharing remains
available as a separate scheme under `nil/crypto3/pubkey/secret_sharing`.

## Scheme And Key Types

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
```

`scheme_type` selects the BLS12-381 curve, the message-signature ciphersuite,
and the basic BLS scheme. The private and public key wrappers expose the types
expected by the generic signing algorithms.

## Sign And Verify

```cpp
int main() {
    const std::string text = "hello world";
    const std::vector<std::uint8_t> message(text.begin(), text.end());

    private_key_type private_key(
        random_element<typename scalar_type::field_type>());
    signature_type signature = sign(message, private_key);

    public_key_type &public_key = private_key;
    assert(verify(message, signature, public_key));
}
```

`sign` and `verify` accept ranges, so the message can be supplied as a byte
container. The private-key wrapper also provides access to its corresponding
public key.

## Build The Maintained Example

From the repository root:

```sh
cmake -S . -B build -DBUILD_EXAMPLES=ON
cmake --build build --target pubkey_bls_example --parallel
./build/libs/pubkey/example/pubkey_bls_example
```

The executable exits silently with status zero when its verification assertion
succeeds.
