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

#ifndef CRYPTO3_MATH_COMPLETE_FACTORIZATION_HPP
#define CRYPTO3_MATH_COMPLETE_FACTORIZATION_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include <nil/crypto3/math/polynomial/equal_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/kaltofen_shoup_distinct_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/square_free_factorization.hpp>

namespace nil::crypto3::math {

    /**
     * Factor a polynomial into monic irreducible factors by composing the three finite-field factorization stages:
     *
     * 1. Square-free factorization separates factors by their multiplicity in the input.
     * 2. Kaltofen-Shoup distinct-degree factorization separates each square-free part into groups whose irreducible
     *    factors all have the same degree.
     * 3. Cantor-Zassenhaus equal-degree factorization splits each group into individual irreducible factors.
     *
     * Thus, a complete result satisfies
     *
     *     input = leading_coefficient * product(factor.polynomial ^ factor.multiplicity).
     *
     * The multiplicity attached to each irreducible factor is inherited from its square-free part. The caller-owned
     * generator is shared by all Cantor-Zassenhaus groups and must return independent uniformly distributed
     * coefficient-field elements.
     *
     * After each irreducible factor is appended to the result, factor_callback may request an early stop. A stopped
     * result includes that factor and has complete set to false. Zero and constant inputs produce no factors and
     * preserve their scalar value as leading_coefficient.
     *
     * @throws std::invalid_argument under the restrictions documented by the component stages, including when the
     * coefficient-field characteristic is not greater than the input degree or is two.
     * @pre input is a nonempty coefficient polynomial.
     */
    template<detail::SupportsDivrem Backend, typename Generator, typename FactorCallback>
        requires detail::PolynomialFactorCallback<FactorCallback, typename Backend::polynomial_type>
    polynomial_factorization_result<typename Backend::polynomial_type>
        complete_factorization(const typename Backend::polynomial_type &input,
                               polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                               Generator &generator, FactorCallback &&factor_callback) {
        using polynomial_type = typename Backend::polynomial_type;
        using result_type = polynomial_factorization_result<polynomial_type>;

        auto square_free_result = square_free_factorization<Backend>(input, arithmetic_context);

        result_type result;
        result.leading_coefficient = square_free_result.leading_coefficient;
        for (auto &square_free_factor : square_free_result.factors) {
            const std::size_t multiplicity = square_free_factor.multiplicity;
            const std::size_t block_size = detail::kaltofen_shoup_block_size(square_free_factor.polynomial.size() - 1);

            std::vector<distinct_degree_factor<polynomial_type>> degree_groups;
            detail::kaltofen_shoup_factor_monic_square_free(degree_groups, std::move(square_free_factor.polynomial),
                                                            block_size, arithmetic_context,
                                                            [](const distinct_degree_factor<polynomial_type> &) {
                                                                return factorization_control::continue_factorization;
                                                            });

            for (auto &degree_group : degree_groups) {
                const factorization_control control = detail::factor_distinct_degree_group<Backend>(
                    std::move(degree_group), arithmetic_context, generator, [&](polynomial_type &&factor) {
                        result.factors.push_back({std::move(factor), multiplicity});
                        return factor_callback(result.factors.back());
                    });
                if (control == factorization_control::stop_factorization) {
                    result.complete = false;
                    return result;
                }
            }
        }
        return result;
    }

    /** Compute the complete irreducible factorization without a staged callback. */
    template<detail::SupportsDivrem Backend, typename Generator>
    polynomial_factorization_result<typename Backend::polynomial_type>
        complete_factorization(const typename Backend::polynomial_type &input,
                               polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                               Generator &generator) {
        using factor_type = polynomial_factor<typename Backend::polynomial_type>;
        return complete_factorization<Backend>(input, arithmetic_context, generator, [](const factor_type &) {
            return factorization_control::continue_factorization;
        });
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_COMPLETE_FACTORIZATION_HPP
