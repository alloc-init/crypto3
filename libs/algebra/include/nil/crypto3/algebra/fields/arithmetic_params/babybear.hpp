//---------------------------------------------------------------------------//
// Copyright (c) 2024 Alexey Kokoshnikov <alexeikokoshnikov@nil.foundation>
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

#ifndef CRYPTO3_ALGEBRA_FIELDS_BABYBEAR_ARITHMETIC_PARAMS_HPP
#define CRYPTO3_ALGEBRA_FIELDS_BABYBEAR_ARITHMETIC_PARAMS_HPP

#include <nil/crypto3/algebra/fields/params.hpp>

#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {

                template<>
                struct arithmetic_params<babybear_base_field> : public params<babybear_base_field> {
                private:
                    typedef params<babybear_base_field> policy_type;

                public:
                    typedef typename policy_type::modular_type modular_type;
                    typedef typename policy_type::integral_type integral_type;

                    constexpr static const std::size_t s = 27;
                    constexpr static integral_type arithmetic_generator = 1u;
                    constexpr static integral_type multiplicative_generator = 31u;
                    constexpr static integral_type root_of_unity = 0x1a427a41u;

                    constexpr static integral_type geometric_generator = 0u;
                };

                constexpr std::size_t const arithmetic_params<babybear_base_field>::s;

                constexpr typename arithmetic_params<babybear_base_field>::integral_type const
                    arithmetic_params<babybear_base_field>::root_of_unity;

                constexpr typename arithmetic_params<babybear_base_field>::integral_type const
                    arithmetic_params<babybear_base_field>::arithmetic_generator;

                constexpr typename arithmetic_params<babybear_base_field>::integral_type const
                    arithmetic_params<babybear_base_field>::geometric_generator;

                constexpr typename arithmetic_params<babybear_base_field>::integral_type const
                    arithmetic_params<babybear_base_field>::multiplicative_generator;

                template<>
                struct arithmetic_params<babybear_fp4> : public params<babybear> {
                    // It's actually 29 but that requires going into the extension field and we don't
                    // want that
                    constexpr static std::size_t s = 27;
                    constexpr static babybear_fp4::value_type multiplicative_generator {
                        {babybear::value_type(8), babybear::value_type(1), babybear::value_type(0),
                         babybear::value_type(0)}};
                    constexpr static integral_type root_of_unity = 0x1a427a41u;

                    constexpr static integral_type arithmetic_generator = 0u;
                    constexpr static integral_type geometric_generator = 0u;
                };
            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ALGEBRA_FIELDS_BABYBEAR_ARITHMETIC_PARAMS_HPP
