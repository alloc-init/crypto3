//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_KDF_DERIVE_HPP
#define CRYPTO3_KDF_DERIVE_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

namespace nil {
    namespace crypto3 {
        // This file is the public algorithm surface for KDFs. Different KDF families have genuinely different input
        // shapes, so future algorithms should add constrained derive<Kdf>(...) overloads instead of forcing every KDF
        // into the PBKDF2 password/salt/iterations/output-length signature.
        /*!
         * @brief Derive key material with a password-based key derivation function.
         *
         * @tparam PasswordBasedKeyDerivationFunction PBKDF policy, for example
         * kdf::pbkdf2<mac::hmac<hashes::sha1>>.
         * @ingroup kdf_algorithms
         */
        template<typename PasswordBasedKeyDerivationFunction, typename PasswordIterator, typename SaltIterator,
                 typename OutputIterator>
        OutputIterator derive(PasswordIterator password_first,
                              PasswordIterator password_last,
                              SaltIterator salt_first,
                              SaltIterator salt_last,
                              std::size_t iterations,
                              std::size_t derived_key_octets,
                              OutputIterator out) {
            return PasswordBasedKeyDerivationFunction::derive(password_first, password_last, salt_first, salt_last,
                                                              iterations, derived_key_octets, out);
        }

        template<typename PasswordBasedKeyDerivationFunction, typename PasswordRange, typename SaltRange,
                 typename OutputIterator>
        OutputIterator derive(const PasswordRange &password,
                              const SaltRange &salt,
                              std::size_t iterations,
                              std::size_t derived_key_octets,
                              OutputIterator out) {
            return PasswordBasedKeyDerivationFunction::derive(password, salt, iterations, derived_key_octets, out);
        }

        template<typename PasswordBasedKeyDerivationFunction, typename PasswordRange, typename SaltRange>
        std::vector<std::uint8_t> derive(const PasswordRange &password,
                                         const SaltRange &salt,
                                         std::size_t iterations,
                                         std::size_t derived_key_octets) {
            return PasswordBasedKeyDerivationFunction::derive(password, salt, iterations, derived_key_octets);
        }
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_KDF_DERIVE_HPP
