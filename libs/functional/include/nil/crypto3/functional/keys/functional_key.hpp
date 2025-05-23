//---------------------------------------------------------------------------//
// Copyright (c) 2024 Mikhail Komarov <nemo@allocin.it>
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

#ifndef CRYPTO3_FUNCTIONAL_PUBLIC_KEY_HPP
#define CRYPTO3_FUNCTIONAL_PUBLIC_KEY_HPP

namespace nil {
    namespace crypto3 {
        namespace functional {
            /*!
             * @brief
             *
             * @ingroup pubkey_algorithms
             *
             * Public key - a key that can be published and used to verify the authenticity
             * of the signed document, as well as to prevent fraud on the part of the certifying
             * person in the form of his refusal to sign the document.
             *
             */
            template<typename SchemeType, typename T = void>
            struct functional_key {
                typedef SchemeType scheme_type;
            };
        } // namespace pubkey
    } // namespace crypto3
} // namespace nil

#endif    // CRYPTO3_PUBKEY_PUBLIC_KEY_HPP
