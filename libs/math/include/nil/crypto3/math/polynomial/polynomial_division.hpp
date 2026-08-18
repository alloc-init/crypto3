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
#include <memory>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/math/polynomial/power_series.hpp>

namespace nil::crypto3::math {

    namespace detail {
        inline bool use_basecase_division(const polynomial_arithmetic::polynomial_context_options &options,
                                          std::size_t divisor_coefficient_count,
                                          std::size_t quotient_coefficient_count) {
            return (options.basecase_divisor_coefficient_cutoff != 0 &&
                    divisor_coefficient_count <= options.basecase_divisor_coefficient_cutoff) ||
                   (options.basecase_quotient_coefficient_cutoff != 0 &&
                    quotient_coefficient_count <= options.basecase_quotient_coefficient_cutoff);
        }
    }    // namespace detail

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
    class polynomial_divisor_context {
    public:
        using backend_type = Backend;
        using polynomial_type = typename backend_type::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        polynomial_divisor_context(const polynomial_type &divisor, std::size_t inverse_precision,
                                   polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context)
            requires std::copy_constructible<polynomial_type> &&
                         requires(polynomial_type &output, const polynomial_type &input, std::size_t coefficient_count,
                                  polynomial_arithmetic::polynomial_context<backend_type> &context) {
                             inverse_series(output, input, coefficient_count, context);
                         }
            : divisor_(divisor), inverse_precision_(inverse_precision) {
            condense(divisor_);
            if (divisor_.size() == 1 && divisor_[0] == value_type {}) {
                throw std::invalid_argument("the zero polynomial cannot be used as a divisor");
            }
            if (inverse_precision_ == 0) {
                throw std::invalid_argument("the divisor inverse precision must be positive");
            }

            polynomial_type reversed_divisor(divisor_);
            reverse(reversed_divisor, reversed_divisor.size());
            inverse_series(reversed_divisor_inverse_, reversed_divisor, inverse_precision_, arithmetic_context);
        }

        const polynomial_type &divisor() const {
            return divisor_;
        }

        std::size_t degree() const {
            return divisor_.size() - 1;
        }

        const polynomial_type &reversed_divisor_inverse() const {
            return reversed_divisor_inverse_;
        }

        std::size_t inverse_precision() const {
            return inverse_precision_;
        }

    private:
        polynomial_type divisor_;
        polynomial_type reversed_divisor_inverse_;
        std::size_t inverse_precision_;
    };

    /**
     * Divide dividend A by the divisor B stored in divisor_context and store the canonical quotient Q and remainder R.
     * The quotient and remainder must be distinct, but either may alias the dividend.
     *
     * For n = degree(A), d = degree(B), and k = n - d + 1, quotient reversal turns division into multiplication:
     *
     *     rev(Q) = rev(A) * rev(B)^-1 mod X^k.
     *
     * After recovering Q, only the first d coefficients of A - Q * B are needed because the remainder has degree less
     * than d. If n < d, the quotient is zero and the dividend is returned unchanged as the remainder.
     *
     * For small divisors or quotients, the arithmetic context selects quadratic long division instead. Its two
     * inclusive coefficient-count cutoffs are independently configurable, and setting either cutoff to zero disables
     * that criterion. The defaults keep very small operations out of the Newton path.
     *
     * @throws std::invalid_argument if quotient and remainder are the same object or the inverse was precomputed to
     *         precision less than k.
     * @pre dividend is a nonempty canonical coefficient polynomial.
     */
    template<polynomial_arithmetic::PolynomialBackend Backend>
        requires detail::MutableNormalizableCoefficientPolynomial<typename Backend::polynomial_type> &&
                 std::default_initializable<typename Backend::polynomial_type> &&
                 std::movable<typename Backend::polynomial_type> &&
                 requires(typename Backend::polynomial_type &quotient, typename Backend::polynomial_type &remainder,
                          const typename Backend::polynomial_type &dividend,
                          const typename Backend::polynomial_type &divisor) {
                     division(quotient, remainder, dividend, divisor);
                 } &&
                 requires(const typename Backend::polynomial_type::value_type &left,
                          const typename Backend::polynomial_type::value_type &right) {
                     { left - right } -> std::convertible_to<typename Backend::polynomial_type::value_type>;
                 }
    void divrem(typename Backend::polynomial_type &quotient, typename Backend::polynomial_type &remainder,
                const typename Backend::polynomial_type &dividend,
                const polynomial_divisor_context<Backend> &divisor_context,
                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if (std::addressof(quotient) == std::addressof(remainder)) {
            throw std::invalid_argument("quotient and remainder must be distinct objects");
        }

        polynomial_type quotient_result;
        polynomial_type remainder_result;
        if (dividend.size() < divisor_context.divisor().size()) {
            quotient_result.resize(1);
            quotient_result[0] = value_type {};
            remainder_result = dividend;
            quotient = std::move(quotient_result);
            remainder = std::move(remainder_result);
            return;
        }

        const std::size_t quotient_size = dividend.size() - divisor_context.divisor().size() + 1;
        const auto &options = arithmetic_context.options();
        // Long division avoids Newton multiplication overhead when either the divisor or quotient is small. It does
        // not use the precomputed reversed-divisor inverse, so this dispatch precedes the inverse-precision check.
        if (detail::use_basecase_division(options, divisor_context.divisor().size(), quotient_size)) {
            division(quotient_result, remainder_result, dividend, divisor_context.divisor());
            quotient = std::move(quotient_result);
            remainder = std::move(remainder_result);
            return;
        }

        if (quotient_size > divisor_context.inverse_precision()) {
            throw std::invalid_argument("the precomputed divisor inverse has insufficient precision");
        }

        polynomial_type reversed_dividend;
        reversed_dividend.resize(quotient_size);
        for (std::size_t i = 0; i < quotient_size; ++i) {
            reversed_dividend[i] = dividend[dividend.size() - 1 - i];
        }

        polynomial_type reversed_quotient;
        arithmetic_context.multiply_low(reversed_quotient, reversed_dividend,
                                        divisor_context.reversed_divisor_inverse(), quotient_size);
        reversed_quotient.resize(quotient_size, value_type {});
        reverse(reversed_quotient, quotient_size);
        condense(reversed_quotient);
        quotient_result = std::move(reversed_quotient);

        const std::size_t divisor_degree = divisor_context.degree();
        if (divisor_degree == 0) {
            remainder_result.resize(1);
            remainder_result[0] = value_type {};
        } else {
            arithmetic_context.multiply_low(remainder_result, quotient_result, divisor_context.divisor(),
                                            divisor_degree);
            // multiply_low returns canonical output; restore the complete prefix before coefficient-wise subtraction.
            remainder_result.resize(divisor_degree, value_type {});
            for (std::size_t i = 0; i < divisor_degree; ++i) {
                remainder_result[i] = dividend[i] - remainder_result[i];
            }
            condense(remainder_result);
        }

        quotient = std::move(quotient_result);
        remainder = std::move(remainder_result);
    }

    namespace detail {
        template<typename Backend>
        concept SupportsDivrem =
            polynomial_arithmetic::PolynomialBackend<Backend> &&
            std::default_initializable<typename Backend::polynomial_type> &&
            requires(typename Backend::polynomial_type &quotient, typename Backend::polynomial_type &remainder,
                     const typename Backend::polynomial_type &dividend,
                     const polynomial_divisor_context<Backend> &divisor_context,
                     polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
                divrem(quotient, remainder, dividend, divisor_context, arithmetic_context);
            };
    }    // namespace detail

    /**
     * Reduce dividend modulo the divisor stored in divisor_context. The result is canonical and may alias dividend.
     * The precomputed inverse must have enough precision for the quotient that divrem computes internally.
     */
    template<detail::SupportsDivrem Backend>
    void remainder(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &dividend,
                   const polynomial_divisor_context<Backend> &divisor_context,
                   polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        typename Backend::polynomial_type quotient;
        divrem(quotient, output, dividend, divisor_context, arithmetic_context);
    }

    /**
     * Divide dividend by the divisor stored in divisor_context and reject a nonzero remainder. The canonical quotient
     * may alias dividend. Output is replaced only after exact divisibility has been verified.
     *
     * @throws std::invalid_argument if dividend is not exactly divisible by the stored divisor or the precomputed
     *         inverse has insufficient precision.
     */
    template<detail::SupportsDivrem Backend>
    void exact_division(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &dividend,
                        const polynomial_divisor_context<Backend> &divisor_context,
                        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        typename Backend::polynomial_type quotient;
        typename Backend::polynomial_type remainder_result;
        divrem(quotient, remainder_result, dividend, divisor_context, arithmetic_context);
        if (!is_zero(remainder_result)) {
            throw std::invalid_argument("polynomial division is not exact");
        }
        output = std::move(quotient);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_DIVISION_HPP
