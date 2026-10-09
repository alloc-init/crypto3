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

#include <nil/crypto3/math/polynomial/arithmetic/power_series.hpp>

namespace nil::crypto3::math {

    namespace detail {
        template<typename Backend>
        concept SupportsPreparedLowProduct =
            requires(polynomial_arithmetic::polynomial_context<Backend> &context,
                     typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                     std::size_t count) {
                { context.prepare_low_product(input, count, count).has_value() } -> std::same_as<bool>;
                {
                    context.try_multiply_low_prepared(output, input, *context.prepare_low_product(input, count, count),
                                                      count)
                } -> std::same_as<bool>;
            };

        // Keep the optional spectra out of backends that only implement the required multiplication interface.
        template<typename Backend>
        struct divisor_prepared_products { };

        template<SupportsPreparedLowProduct Backend>
        struct divisor_prepared_products<Backend> {
            using preparation_type =
                decltype(std::declval<polynomial_arithmetic::polynomial_context<Backend> &>().prepare_low_product(
                    std::declval<const typename Backend::polynomial_type &>(), std::size_t {}, std::size_t {}));
            preparation_type inverse;
            preparation_type divisor;
        };

        template<typename Backend>
        concept SupportsPreparedCyclicRemainder =
            requires(polynomial_arithmetic::polynomial_context<Backend> &context,
                     typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                     std::size_t count) {
                { context.prepare_cyclic_remainder(input, count).has_value() } -> std::same_as<bool>;
                {
                    context.try_cyclic_remainder(output, input, input, *context.prepare_cyclic_remainder(input, count))
                } -> std::same_as<bool>;
            };

        template<typename Backend>
        struct divisor_prepared_cyclic_remainder { };

        template<SupportsPreparedCyclicRemainder Backend>
        struct divisor_prepared_cyclic_remainder<Backend> {
            using preparation_type =
                decltype(std::declval<polynomial_arithmetic::polynomial_context<Backend> &>().prepare_cyclic_remainder(
                    std::declval<const typename Backend::polynomial_type &>(),
                    std::size_t {}));
            preparation_type value;
        };

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
     * When supported by the backend, also prepare the two fixed low-product operands for a quotient with
     * inverse_precision coefficients. These preparations belong to this divisor snapshot and never grow or change
     * during division. Each product falls back to ordinary multiplication when its precision or transform length
     * is incompatible with its preparation.
     * No preparations are needed when that quotient size selects long division, or when B is constant.
     * The modulus preparation uses a shorter cyclic FFT when supported and beneficial; otherwise it uses the
     * original low-product FFT. At most one of these modulus spectra is retained, alongside the inverse spectrum.
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

        template<polynomial_arithmetic::PolynomialObserver Observer>
        polynomial_divisor_context(
            const polynomial_type &divisor, std::size_t inverse_precision,
            polynomial_arithmetic::polynomial_context<backend_type, Observer> &arithmetic_context)
            requires std::copy_constructible<polynomial_type> &&
                         requires(polynomial_type &output, const polynomial_type &input, std::size_t coefficient_count,
                                  polynomial_arithmetic::polynomial_context<backend_type, Observer> &context) {
                             inverse_series(output, input, coefficient_count, context);
                         }
            : divisor_(divisor), inverse_precision_(inverse_precision) {
            using stage = polynomial_arithmetic::polynomial_stage;
            using metric = polynomial_arithmetic::polynomial_metric;
            using metadata = polynomial_arithmetic::polynomial_metadata;
            // Member initialization has already copied the owned divisor; this scope covers its precomputation.
            auto scope = arithmetic_context.template observe<stage::divisor_preparation>([&]() noexcept {
                return metadata {
                    {{metric::input_coefficients, divisor_.size()}, {metric::inverse_precision, inverse_precision_}}};
            });
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
            if constexpr (detail::SupportsPreparedLowProduct<backend_type>) {
                if (degree() != 0 &&
                    !detail::use_basecase_division(arithmetic_context.options(), divisor_.size(), inverse_precision_)) {
                    prepared_.inverse = arithmetic_context.prepare_low_product(reversed_divisor_inverse_,
                                                                               inverse_precision_, inverse_precision_);
                    if constexpr (detail::SupportsPreparedCyclicRemainder<backend_type>) {
                        cyclic_.value = arithmetic_context.prepare_cyclic_remainder(divisor_, inverse_precision_);
                        if (cyclic_.value) {
                            return;
                        }
                    }
                    prepared_.divisor = arithmetic_context.prepare_low_product(divisor_, inverse_precision_, degree());
                }
            }
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

        // Multiply by the stored reversed inverse, retaining exactly coefficient_count low coefficients. The
        // preparation must match this precision and the selected FFT length; otherwise use the ordinary product.
        template<polynomial_arithmetic::PolynomialObserver Observer>
        void multiply_low_reversed_inverse(
            polynomial_type &output, const polynomial_type &input, std::size_t coefficient_count,
            polynomial_arithmetic::polynomial_context<backend_type, Observer> &context) const {
            if constexpr (detail::SupportsPreparedLowProduct<backend_type>) {
                if (prepared_.inverse &&
                    context.try_multiply_low_prepared(output, input, *prepared_.inverse, coefficient_count)) {
                    return;
                }
            }
            // This is the existing fallback, either without preparation or after a declined prepared attempt.
            auto fallback =
                context.template observe<polynomial_arithmetic::polynomial_stage::low_product_fallback>([&]() noexcept {
                    using metric = polynomial_arithmetic::polynomial_metric;
                    polynomial_arithmetic::polynomial_metadata result {
                        {{metric::low_coefficient_count, coefficient_count}, {metric::prepared_available, 0}}};
                    if constexpr (detail::SupportsPreparedLowProduct<backend_type>) {
                        result[1].value = bool(prepared_.inverse);
                    }
                    return result;
                });
            context.multiply_low(output, input, reversed_divisor_inverse_, coefficient_count);
        }

        // Only degree(B) low coefficients of Q*B are needed for remainder reconstruction. The leading coefficient
        // of B cannot affect this prefix. Keep the same low-product operation whether or not preparation succeeds.
        template<polynomial_arithmetic::PolynomialObserver Observer>
        void multiply_low_divisor(polynomial_type &output, const polynomial_type &input,
                                  polynomial_arithmetic::polynomial_context<backend_type, Observer> &context) const {
            if constexpr (detail::SupportsPreparedLowProduct<backend_type>) {
                if (prepared_.divisor &&
                    context.try_multiply_low_prepared(output, input, *prepared_.divisor, degree())) {
                    return;
                }
            }
            // This is the existing fallback, either without preparation or after a declined prepared attempt.
            auto fallback =
                context.template observe<polynomial_arithmetic::polynomial_stage::low_product_fallback>([&]() noexcept {
                    using metric = polynomial_arithmetic::polynomial_metric;
                    polynomial_arithmetic::polynomial_metadata result {
                        {{metric::low_coefficient_count, degree()}, {metric::prepared_available, 0}}};
                    if constexpr (detail::SupportsPreparedLowProduct<backend_type>) {
                        result[1].value = bool(prepared_.divisor);
                    }
                    return result;
                });
            context.multiply_low(output, input, divisor_, degree());
        }

        // With exact quotient Q, cyclic reconstruction cancels the wrapped high coefficients of V and Q*B.
        // The backend checks its bounded degree range and current transform lengths before changing output.
        template<polynomial_arithmetic::PolynomialObserver Observer>
        bool try_reconstruct_remainder(
            polynomial_type &output, const polynomial_type &dividend, const polynomial_type &quotient,
            polynomial_arithmetic::polynomial_context<backend_type, Observer> &context) const {
            if constexpr (detail::SupportsPreparedCyclicRemainder<backend_type>) {
                if (cyclic_.value) {
                    return context.try_cyclic_remainder(output, dividend, quotient, *cyclic_.value);
                }
            }
            return false;
        }

    private:
        polynomial_type divisor_;
        polynomial_type reversed_divisor_inverse_;
        std::size_t inverse_precision_;
        [[no_unique_address]] detail::divisor_prepared_products<backend_type> prepared_;
        [[no_unique_address]] detail::divisor_prepared_cyclic_remainder<backend_type> cyclic_;
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
     * For n <= 2*d - 2, a backend may reconstruct R using a shorter cyclic convolution instead: folding A and
     * Q*B modulo X^M - 1 with M >= d cancels their wrapped high coefficients and leaves exactly R.
     *
     * For small divisors or quotients, the arithmetic context selects quadratic long division instead. Its two
     * inclusive coefficient-count cutoffs are independently configurable, and setting either cutoff to zero disables
     * that criterion. The defaults keep very small operations out of the Newton path.
     *
     * @throws std::invalid_argument if quotient and remainder are the same object or the inverse was precomputed to
     *         precision less than k.
     * @pre dividend is a nonempty canonical coefficient polynomial.
     */
    template<polynomial_arithmetic::PolynomialBackend Backend, polynomial_arithmetic::PolynomialObserver Observer>
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
                polynomial_arithmetic::polynomial_context<Backend, Observer> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        using stage = polynomial_arithmetic::polynomial_stage;
        using metric = polynomial_arithmetic::polynomial_metric;
        using metadata = polynomial_arithmetic::polynomial_metadata;
        using path = polynomial_arithmetic::polynomial_arithmetic_path;
        auto scope = arithmetic_context.template observe<stage::division>([&]() noexcept {
            return metadata {
                {{metric::input_coefficients, dividend.size()}, {metric::modulus_degree, divisor_context.degree()}}};
        });

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
            scope.set_result_metadata([&]() noexcept {
                return metadata {{{metric::quotient_coefficients, quotient.size()},
                                  {metric::remainder_coefficients, remainder.size()},
                                  {metric::arithmetic_path, static_cast<std::size_t>(path::trivial)}}};
            });
            return;
        }

        const std::size_t quotient_size = dividend.size() - divisor_context.divisor().size() + 1;
        const auto &options = arithmetic_context.options();
        // Long division avoids Newton multiplication overhead when either the divisor or quotient is small. It does
        // not use the precomputed reversed-divisor inverse, so this dispatch precedes the inverse-precision check.
        if (detail::use_basecase_division(options, divisor_context.divisor().size(), quotient_size)) {
            auto basecase = arithmetic_context.template observe<stage::basecase_division>();
            division(quotient_result, remainder_result, dividend, divisor_context.divisor());
            quotient = std::move(quotient_result);
            remainder = std::move(remainder_result);
            scope.set_result_metadata([&]() noexcept {
                return metadata {{{metric::quotient_coefficients, quotient.size()},
                                  {metric::remainder_coefficients, remainder.size()},
                                  {metric::arithmetic_path, static_cast<std::size_t>(path::basecase)}}};
            });
            return;
        }

        if (quotient_size > divisor_context.inverse_precision()) {
            throw std::invalid_argument("the precomputed divisor inverse has insufficient precision");
        }

        // Each estimation reuses the coefficient-space inverse; its construction is a separate scope.
        auto estimate = arithmetic_context.template observe<stage::quotient_estimation>([&]() noexcept {
            return metadata {{{metric::quotient_coefficients, quotient_size},
                              {metric::inverse_precision, divisor_context.inverse_precision()}}};
        });
        polynomial_type reversed_dividend;
        reversed_dividend.resize(quotient_size);
        for (std::size_t i = 0; i < quotient_size; ++i) {
            reversed_dividend[i] = dividend[dividend.size() - 1 - i];
        }

        polynomial_type reversed_quotient;
        divisor_context.multiply_low_reversed_inverse(reversed_quotient, reversed_dividend, quotient_size,
                                                      arithmetic_context);
        reversed_quotient.resize(quotient_size, value_type {});
        reverse(reversed_quotient, quotient_size);
        condense(reversed_quotient);
        quotient_result = std::move(reversed_quotient);
        estimate.finish();

        {
            auto reconstruction = arithmetic_context.template observe<stage::remainder_reconstruction>([&]() noexcept {
                return metadata {{{metric::input_coefficients, dividend.size()},
                                  {metric::quotient_coefficients, quotient_result.size()},
                                  {metric::modulus_degree, divisor_context.degree()}}};
            });
            const std::size_t divisor_degree = divisor_context.degree();
            if (divisor_degree == 0) {
                remainder_result.resize(1);
                remainder_result[0] = value_type {};
                reconstruction.set_result_metadata([]() noexcept {
                    return metadata {{{metric::arithmetic_path, static_cast<std::size_t>(path::trivial)}}};
                });
            } else if (!divisor_context.try_reconstruct_remainder(remainder_result, dividend, quotient_result,
                                                                  arithmetic_context)) {
                divisor_context.multiply_low_divisor(remainder_result, quotient_result, arithmetic_context);
                // multiply_low returns canonical output; restore the complete prefix before coefficient-wise
                // subtraction.
                remainder_result.resize(divisor_degree, value_type {});
                for (std::size_t i = 0; i < divisor_degree; ++i) {
                    remainder_result[i] = dividend[i] - remainder_result[i];
                }
                condense(remainder_result);
                reconstruction.set_result_metadata([&]() noexcept {
                    return metadata {{{metric::arithmetic_path, static_cast<std::size_t>(path::low_product)},
                                      {metric::output_coefficients, remainder_result.size()}}};
                });
            } else {
                reconstruction.set_result_metadata([&]() noexcept {
                    return metadata {{{metric::arithmetic_path, static_cast<std::size_t>(path::cyclic)},
                                      {metric::output_coefficients, remainder_result.size()}}};
                });
            }
        }

        quotient = std::move(quotient_result);
        remainder = std::move(remainder_result);
        scope.set_result_metadata([&]() noexcept {
            return metadata {{{metric::quotient_coefficients, quotient.size()},
                              {metric::remainder_coefficients, remainder.size()},
                              {metric::arithmetic_path, static_cast<std::size_t>(path::reciprocal)}}};
        });
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
    template<detail::SupportsDivrem Backend, polynomial_arithmetic::PolynomialObserver Observer>
    void remainder(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &dividend,
                   const polynomial_divisor_context<Backend> &divisor_context,
                   polynomial_arithmetic::polynomial_context<Backend, Observer> &arithmetic_context) {
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
    template<detail::SupportsDivrem Backend, polynomial_arithmetic::PolynomialObserver Observer>
    void exact_division(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &dividend,
                        const polynomial_divisor_context<Backend> &divisor_context,
                        polynomial_arithmetic::polynomial_context<Backend, Observer> &arithmetic_context) {
        using stage = polynomial_arithmetic::polynomial_stage;
        using metric = polynomial_arithmetic::polynomial_metric;
        using metadata = polynomial_arithmetic::polynomial_metadata;
        auto scope = arithmetic_context.template observe<stage::exact_division>([&]() noexcept {
            return metadata {
                {{metric::input_coefficients, dividend.size()}, {metric::modulus_degree, divisor_context.degree()}}};
        });
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
