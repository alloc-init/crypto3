//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021 Nikita Kaskov <nbering@nil.foundation>
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

#ifndef CRYPTO3_MATH_POLYNOMIAL_SHIFT_HPP
#define CRYPTO3_MATH_POLYNOMIAL_SHIFT_HPP

#include <cassert>
#include <concepts>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>

#include <nil/crypto3/algebra/type_traits.hpp>
#include <nil/crypto3/math/polynomial/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>
#include <nil/crypto3/math/polynomial/polynomial_dfs.hpp>

namespace nil {
    namespace crypto3 {
        namespace math {
            /**
             * Scale the polynomial argument: return
             * g(X) = f(x * X) = sum_i a_i * x^i * X^i for f(X) = sum_i a_i * X^i.
             * This changes coefficient values without moving them to different degrees.
             */
            template<typename FieldValueType>
                requires algebra::is_field_element<FieldValueType>::value
            static inline polynomial<FieldValueType> polynomial_shift(const polynomial<FieldValueType> &f,
                                                                      const FieldValueType &x) {
                polynomial<FieldValueType> f_shifted(f);
                FieldValueType x_power = x;
                for (std::size_t i = 1; i < f.size(); i++) {
                    f_shifted[i] *= x_power;
                    if (i + 1 < f.size()) {
                        x_power *= x;
                    }
                }

                return f_shifted;
            }

            /**
             * Shift coefficients toward higher degrees: compute g(X) = X^shift * f(X). With ascending coefficient
             * storage, this prepends shift zero coefficients. Store the canonical result in output, which may alias
             * input.
             */
            template<CoefficientPolynomial Polynomial>
                requires detail::MutablePolynomialCoefficientRange<Polynomial> &&
                         std::default_initializable<typename Polynomial::value_type> &&
                         std::equality_comparable<typename Polynomial::value_type>
            void shift_left(Polynomial &output, const Polynomial &input, std::size_t shift) {
                using value_type = typename Polynomial::value_type;

                if (input.size() == 1 && input[0] == value_type {}) {
                    output.resize(1);
                    output[0] = value_type {};
                    return;
                }

                const std::size_t input_size = input.size();
                if (shift > std::numeric_limits<std::size_t>::max() - input_size) {
                    throw std::length_error("polynomial shift exceeds the maximum coefficient count");
                }
                output.resize(input_size + shift, value_type {});
                if (std::addressof(output) == std::addressof(input)) {
                    for (std::size_t i = input_size; i > 0; --i) {
                        output[i - 1 + shift] = output[i - 1];
                    }
                } else {
                    for (std::size_t i = 0; i < input_size; ++i) {
                        output[i + shift] = input[i];
                    }
                }
                for (std::size_t i = 0; i < shift; ++i) {
                    output[i] = value_type {};
                }
                condense(output);
            }

            /**
             * Shift coefficients toward lower degrees: for f(X) = sum_i a_i * X^i, compute
             * g(X) = sum_{i >= shift} a_i * X^(i - shift). With ascending coefficient storage, this discards the
             * first shift coefficients. Store the canonical result in output, which may alias input.
             */
            template<CoefficientPolynomial Polynomial>
                requires detail::MutablePolynomialCoefficientRange<Polynomial> &&
                         std::default_initializable<typename Polynomial::value_type> &&
                         std::equality_comparable<typename Polynomial::value_type>
            void shift_right(Polynomial &output, const Polynomial &input, std::size_t shift) {
                using value_type = typename Polynomial::value_type;

                if (shift >= input.size()) {
                    output.resize(1);
                    output[0] = value_type {};
                    return;
                }

                const std::size_t result_size = input.size() - shift;
                if (std::addressof(output) == std::addressof(input)) {
                    for (std::size_t i = 0; i < result_size; ++i) {
                        output[i] = output[i + shift];
                    }
                    output.resize(result_size);
                } else {
                    output.resize(result_size);
                    for (std::size_t i = 0; i < result_size; ++i) {
                        output[i] = input[i + shift];
                    }
                }
                condense(output);
            }

            /**
             * Scale the argument of the represented polynomial by a root of unity: return the DFS representation
             * of g(X) = f(omega^shift * X), where omega is the root for the logical domain_size. This rotates the
             * evaluations cyclically, including when f is stored over a larger extended domain.
             *
             * @pre f is nonempty.
             * @pre domain_size is zero, selecting f.size(), or is a positive divisor of f.size().
             */
            template<typename FieldValueType>
                requires algebra::is_field_element<FieldValueType>::value
            static inline polynomial_dfs<FieldValueType> polynomial_shift(const polynomial_dfs<FieldValueType> &f,
                                                                          const int shift,
                                                                          std::size_t domain_size = 0) {
                if (domain_size == 0) {
                    domain_size = f.size();
                }

                const std::size_t extended_domain_size = f.size();

                assert((extended_domain_size % domain_size) == 0);

                const std::size_t domain_scale = extended_domain_size / domain_size;
                std::size_t normalized_shift;
                if (shift >= 0) {
                    normalized_shift = static_cast<std::size_t>(shift) % domain_size;
                } else {
                    const std::size_t magnitude = static_cast<std::size_t>(-(static_cast<std::int64_t>(shift)));
                    const std::size_t remainder = magnitude % domain_size;
                    normalized_shift = remainder == 0 ? 0 : domain_size - remainder;
                }

                polynomial_dfs<FieldValueType> f_shifted(f.degree(), extended_domain_size);

                for (std::size_t index = 0; index < extended_domain_size; ++index) {
                    f_shifted[index] = f[(index + domain_scale * normalized_shift) % extended_domain_size];
                }

                return f_shifted;
            }
        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_PLONK_REDSHIFT_POLYNOMIAL_SHIFT_HPP
