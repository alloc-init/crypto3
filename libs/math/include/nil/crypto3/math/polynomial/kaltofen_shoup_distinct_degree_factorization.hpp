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

#ifndef CRYPTO3_MATH_KALTOFEN_SHOUP_DISTINCT_DEGREE_FACTORIZATION_HPP
#define CRYPTO3_MATH_KALTOFEN_SHOUP_DISTINCT_DEGREE_FACTORIZATION_HPP

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/math/polynomial/polynomial_frobenius.hpp>

namespace nil::crypto3::math::detail {

    /**
     * Frobenius precomputation for one Kaltofen-Shoup degree block. Let the coefficient field contain Q elements, let
     * B be the modulus polynomial stored in frobenius_context, and let b be the degree-block size. The precomputation
     * stores the baby steps
     *
     *     h[i] = X^(Q^i) mod B,  0 <= i <= b,
     *
     * and a modular-composition precomputation for h[b]. Composing a reduced polynomial A with h[b] applies b
     * Frobenius iterations at once:
     *
     *     A(h[b]) mod B = A(X)^(Q^b) mod B.
     *
     * Reusing that composition advances consecutive giant steps with one modular composition per block.
     *
     * Constructing the baby-step table performs b Frobenius maps and stores b + 1 residues, using
     * O(b * degree(B)) field elements. Construction also builds the existing modular-composition precomputation for
     * h[b].
     *
     * @throws std::invalid_argument if block_size is zero.
     */
    template<SupportsDivrem Backend>
    class kaltofen_shoup_frobenius_precomputation {
    public:
        using backend_type = Backend;
        using polynomial_type = typename backend_type::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        kaltofen_shoup_frobenius_precomputation(
            std::size_t block_size, const polynomial_frobenius_context<backend_type> &frobenius_context,
            polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context) :
            block_size_(block_size), baby_steps_(make_baby_steps(block_size_, frobenius_context, arithmetic_context)),
            block_frobenius_precomputation_(baby_steps_.back(),
                                            // A residue modulo degree-d polynomial B has at most d coefficients. The
                                            // composition precomputation requires a positive limit for constant B.
                                            std::max<std::size_t>(1, frobenius_context.divisor_context().degree()),
                                            frobenius_context.divisor_context(), arithmetic_context) {
        }

        std::size_t block_size() const {
            return block_size_;
        }

        /** Return h[index]. @pre index <= block_size(). */
        const polynomial_type &baby_step(std::size_t index) const {
            return baby_steps_[index];
        }

        /**
         * Apply b Frobenius iterations with one modular composition. Output may alias input.
         *
         * @pre input is reduced modulo B.
         * @pre frobenius_context represents the same B used to construct this precomputation.
         */
        void apply_block_frobenius(polynomial_type &output, const polynomial_type &input,
                                   const polynomial_frobenius_context<backend_type> &frobenius_context,
                                   polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context) const {
            compose_mod(output, input, block_frobenius_precomputation_, frobenius_context.divisor_context(),
                        arithmetic_context);
        }

    private:
        static std::vector<polynomial_type>
            make_baby_steps(std::size_t block_size, const polynomial_frobenius_context<backend_type> &frobenius_context,
                            polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context) {
            if (block_size == 0) {
                throw std::invalid_argument("the Kaltofen-Shoup block size must be positive");
            }

            std::vector<polynomial_type> baby_steps;
            baby_steps.reserve(block_size + 1);
            polynomial_type reduced_x = {value_type {}, value_type::one()};
            const auto &divisor_context = frobenius_context.divisor_context();
            if (divisor_context.degree() == 0) {
                reduced_x.assign(1, value_type {});
            } else if (reduced_x.size() >= divisor_context.divisor().size()) {
                remainder(reduced_x, reduced_x, divisor_context, arithmetic_context);
            }
            baby_steps.emplace_back(std::move(reduced_x));
            for (std::size_t index = 1; index <= block_size; ++index) {
                polynomial_type next_step;
                frobenius_map(next_step, baby_steps.back(), frobenius_context, arithmetic_context);
                baby_steps.emplace_back(std::move(next_step));
            }
            return baby_steps;
        }

        std::size_t block_size_;
        std::vector<polynomial_type> baby_steps_;
        polynomial_composition_precomputation<backend_type> block_frobenius_precomputation_;
    };

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_KALTOFEN_SHOUP_DISTINCT_DEGREE_FACTORIZATION_HPP
