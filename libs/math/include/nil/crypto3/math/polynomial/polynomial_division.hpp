//---------------------------------------------------------------------------//
// Copyright (c) 2026
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

#ifndef CRYPTO3_MATH_POLYNOMIAL_DIVISION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_DIVISION_HPP

#include <concepts>
#include <cstddef>
#include <stdexcept>

#include <nil/crypto3/math/polynomial/power_series.hpp>

namespace nil::crypto3::math {

    /**
     * Immutable precomputation for repeated division and reduction by one polynomial B. If d = degree(B), define
     * rev(B) = X^d * B(X^-1), which reverses B's d + 1 coefficients. The context stores B in canonical form and
     * rev(B)^-1 modulo X^inverse_precision. A later division may reuse this inverse when its quotient has at most
     * inverse_precision coefficients.
     *
     * @pre inverse_precision is positive.
     */
    template<polynomial_arithmetic::PolynomialBackend Backend>
        requires detail::MutableNormalizableCoefficientPolynomial<typename Backend::polynomial_type> &&
                 std::default_initializable<typename Backend::polynomial_type>
    class polynomial_modulus_context {
    public:
        using backend_type = Backend;
        using polynomial_type = typename backend_type::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        polynomial_modulus_context(const polynomial_type &modulus, std::size_t inverse_precision,
                                   polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context)
            requires std::copy_constructible<polynomial_type> &&
                         requires(polynomial_type &output, const polynomial_type &input, std::size_t coefficient_count,
                                  polynomial_arithmetic::polynomial_context<backend_type> &context) {
                             inverse_series(output, input, coefficient_count, context);
                         }
            : modulus_(modulus), inverse_precision_(inverse_precision) {
            condense(modulus_);
            if (modulus_.size() == 1 && modulus_[0] == value_type {}) {
                throw std::invalid_argument("the zero polynomial cannot be used as a modulus");
            }
            if (inverse_precision_ == 0) {
                throw std::invalid_argument("the modulus inverse precision must be positive");
            }

            polynomial_type reversed_modulus(modulus_);
            reverse(reversed_modulus, reversed_modulus.size());
            inverse_series(reversed_modulus_inverse_, reversed_modulus, inverse_precision_, arithmetic_context);
        }

        const polynomial_type &modulus() const {
            return modulus_;
        }

        std::size_t degree() const {
            return modulus_.size() - 1;
        }

        const polynomial_type &reversed_modulus_inverse() const {
            return reversed_modulus_inverse_;
        }

        std::size_t inverse_precision() const {
            return inverse_precision_;
        }

    private:
        polynomial_type modulus_;
        polynomial_type reversed_modulus_inverse_;
        std::size_t inverse_precision_;
    };

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_DIVISION_HPP
