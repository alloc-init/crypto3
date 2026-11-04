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

#ifndef CRYPTO3_MATH_POLYNOMIAL_GCD_HPP
#define CRYPTO3_MATH_POLYNOMIAL_GCD_HPP

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <optional>
#include <utility>

#include <nil/crypto3/math/polynomial/arithmetic/half_gcd.hpp>

namespace nil::crypto3::math {

    namespace detail {

        template<typename Polynomial>
        struct normal_euclidean_chain_certificate {
            explicit normal_euclidean_chain_certificate(const Polynomial &first, const Polynomial &second) :
                expected_divisions(first.size() > 1 ? first.size() - 1 : 0),
                normal(first.size() > 1 && second.size() + 1 == first.size()) {
            }

            void record_quotient(std::size_t coefficients) noexcept {
                ++divisions;
                normal = normal && coefficients == 2;
            }

            void finish(const Polynomial &last_nonzero_remainder) {
                // Counting the complete chain is necessary: the last truncated high-part remainder can lose terms
                // that reconstruction restores. A local high-part degree drop alone is not a global certificate.
                if (normal && divisions == expected_divisions && last_nonzero_remainder.size() == 1 &&
                    !is_zero(last_nonzero_remainder)) {
                    resultant_square_class = last_nonzero_remainder[0];
                }
            }

            std::size_t expected_divisions;
            std::size_t divisions = 0;
            bool normal;
            std::optional<typename Polynomial::value_type> resultant_square_class;
        };

        template<SupportsDivrem Backend, typename ChainCertificate>
        void gcd_impl(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &left,
                      const typename Backend::polynomial_type &right,
                      polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                      ChainCertificate &certificate) {
            using polynomial_type = typename Backend::polynomial_type;

            polynomial_type previous_remainder(left);
            polynomial_type current_remainder(right);
            if (is_zero(previous_remainder)) {
                certificate.finish(current_remainder);
                if (is_zero(current_remainder)) {
                    output = std::move(current_remainder);
                } else {
                    make_monic(output, current_remainder);
                }
                return;
            }
            if (is_zero(current_remainder)) {
                certificate.finish(previous_remainder);
                make_monic(output, previous_remainder);
                return;
            }

            if (previous_remainder.size() < current_remainder.size()) {
                std::swap(previous_remainder, current_remainder);
            }

            polynomial_type quotient;
            polynomial_type next_remainder;
            gcd_divrem_step(quotient, next_remainder, previous_remainder, current_remainder, arithmetic_context,
                            certificate);
            previous_remainder = std::move(current_remainder);
            current_remainder = std::move(next_remainder);

            const auto &options = arithmetic_context.options();
            while (!is_zero(current_remainder)) {
                if (options.gcd_half_gcd_cutoff == 0 || current_remainder.size() < options.gcd_half_gcd_cutoff) {
                    gcd_divrem_step(quotient, next_remainder, previous_remainder, current_remainder, arithmetic_context,
                                    certificate);
                    previous_remainder = std::move(current_remainder);
                    current_remainder = std::move(next_remainder);
                    continue;
                }

                polynomial_type reduced_first;
                polynomial_type reduced_second;
                half_gcd_reduce(reduced_first, reduced_second, previous_remainder, current_remainder,
                                options.half_gcd_basecase_cutoff, arithmetic_context, certificate);
                previous_remainder = std::move(reduced_first);
                current_remainder = std::move(reduced_second);
                if (is_zero(current_remainder)) {
                    break;
                }

                gcd_divrem_step(quotient, next_remainder, previous_remainder, current_remainder, arithmetic_context,
                                certificate);
                previous_remainder = std::move(current_remainder);
                current_remainder = std::move(next_remainder);
            }

            certificate.finish(previous_remainder);
            make_monic(output, previous_remainder);
        }

    }    // namespace detail

    /**
     * Compute the monic greatest common divisor of two canonical coefficient polynomials. Small inputs use the
     * Euclidean algorithm. When enabled by the configured cutoff, larger inputs use recursive half-GCD reduction.
     * Polynomial products and division steps use the configured arithmetic backend. Output may alias either input.
     * The GCD of two zero polynomials is the zero polynomial.
     *
     * The Euclidean algorithm divides the larger polynomial A by the smaller polynomial B and replaces the pair
     * (A, B) with (B, A mod B). The remainder has lower degree than B, so repeating this step eventually leaves a zero
     * remainder. This replacement preserves all common divisors: a divisor of A and B also divides A - q * B, the
     * remainder, while a divisor of B and the remainder also divides q * B + remainder, which is A. At termination,
     * the common divisors are therefore exactly the divisors of the last nonzero remainder, making its monic form the
     * GCD. Half-GCD batches several of these same quotient steps into one recursively constructed polynomial-matrix
     * transformation.
     */
    template<detail::SupportsDivrem Backend>
        requires std::copy_constructible<typename Backend::polynomial_type> &&
                 std::movable<typename Backend::polynomial_type> &&
                 requires(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input) {
                     make_monic(output, input);
                 }
    void gcd(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &left,
             const typename Backend::polynomial_type &right,
             polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        detail::no_euclidean_chain_certificate certificate;
        detail::gcd_impl(output, left, right, arithmetic_context, certificate);
    }

    /**
     * Compute the same monic GCD as gcd, optionally certifying a representative of resultant(left, right) modulo
     * nonzero squares. A certificate is available only for a complete normal Euclidean chain: the inputs have
     * degrees n and n - 1, all n quotient steps are linear, and the final nonzero remainder is a constant c.
     * The returned representative is c, captured before monic normalization. Other chains return std::nullopt;
     * absence is not evidence about the resultant's square class or the input's factorization.
     *
     * For consecutive remainder degrees d, d - 1, d - 2, the resultant recurrence contributes the leading
     * coefficient of the middle remainder to power two. Its sign is (-1)^(d*(d - 1)) = 1. At the end,
     * resultant(a*X + b, c) = c, so resultant(left, right)/c is a nonzero square. Counting all quotient steps,
     * including recursive high-part steps, checks the complete chain without storing the remainders.
     *
     * For monic square-free F of degree n over an odd finite field, factor-count parity additionally requires
     * the discriminant sign: (-1)^r = (-1)^n * chi((-1)^(n*(n - 1)/2) * c), using left = F, right = F'.
     * The sign inside chi is a square over BN254 Fp12. This function does not perform the character test and does
     * not infer irreducibility. It has the same canonical-input and aliasing requirements as gcd.
     */
    template<detail::SupportsDivrem Backend>
        requires std::copy_constructible<typename Backend::polynomial_type> &&
                 std::movable<typename Backend::polynomial_type> &&
                 requires(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input) {
                     make_monic(output, input);
                 }
    std::optional<typename Backend::polynomial_type::value_type>
        gcd_with_normal_chain_certificate(typename Backend::polynomial_type &output,
                                          const typename Backend::polynomial_type &left,
                                          const typename Backend::polynomial_type &right,
                                          polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        detail::normal_euclidean_chain_certificate<typename Backend::polynomial_type> certificate(left, right);
        detail::gcd_impl(output, left, right, arithmetic_context, certificate);
        return certificate.resultant_square_class;
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_GCD_HPP
