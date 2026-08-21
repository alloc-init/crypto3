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
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/math/detail/integer_sqrt.hpp>
#include <nil/crypto3/math/polynomial/distinct_degree_factorization.hpp>
#include <nil/crypto3/math/polynomial/gcd.hpp>
#include <nil/crypto3/math/polynomial/polynomial_factorization.hpp>
#include <nil/crypto3/math/polynomial/polynomial_frobenius.hpp>

namespace nil::crypto3::math::detail {

    /**
     * Select the Kaltofen-Shoup degree-block size for a degree-n polynomial. Distinct-degree factorization explicitly
     * examines degrees only through n / 2; any factor left afterward is irreducible. Choosing
     *
     *     b = ceil(sqrt(n / 2))
     *
     * balances the b baby Frobenius steps against approximately n / (2b) giant-step blocks: the sum
     * b + n / (2b) is minimized at b = sqrt(n / 2). The computation rounds n / 2 upward before applying the integer
     * ceiling square root, which gives the same integer b without using floating-point arithmetic. Constant and
     * linear inputs use the minimum valid block size of one.
     */
    constexpr std::size_t kaltofen_shoup_block_size(std::size_t polynomial_degree) {
        if (polynomial_degree <= 1) {
            return 1;
        }
        const std::size_t half_degree_rounded_up = polynomial_degree / 2 + polynomial_degree % 2;
        return ceil_sqrt(half_degree_rounded_up);
    }

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

    /**
     * Extract one coarse Kaltofen-Shoup degree block. Let the coefficient field contain Q elements, let B be the
     * modulus in frobenius_context, and let precomputation store h[i] = X^(Q^i) mod B for a block of size b. For
     * zero-based block j, the giant-step exponent is e = (j + 1) * b. The number c = degree_count satisfies
     * 1 <= c <= b and allows the final block to process fewer than b degrees. Let R denote the polynomial passed as
     * remaining. If giant_step = X^(Q^e) mod B, form
     *
     *     interval_product = product(giant_step - h[b - 1 - i]) mod B,  0 <= i < c.
     *
     * These differences correspond to degrees e - b + 1 through e - b + c. Their product permits one GCD with R to
     * extract every irreducible factor in that interval. A difference may also contain factors whose degrees divide
     * the target degree, so R must already have all lower-degree factors removed.
     *
     * For example, with b = 2, the first giant step is X^(Q^2): subtracting h[1] targets degree one and subtracting
     * h[0] targets degree two. The next giant step is X^(Q^4): the same subtractions target degrees three and four.
     *
     * @throws std::invalid_argument if degree_count is zero or greater than the block size.
     * @pre remaining is square-free, divides B, and has all earlier degree factors removed.
     * @pre giant_step is reduced modulo B.
     * @pre frobenius_context represents the same B used to construct precomputation.
     */
    template<SupportsDivrem Backend>
    void kaltofen_shoup_coarse_block_factor(typename Backend::polynomial_type &output,
                                            const typename Backend::polynomial_type &remaining,
                                            const typename Backend::polynomial_type &giant_step,
                                            std::size_t degree_count,
                                            const kaltofen_shoup_frobenius_precomputation<Backend> &precomputation,
                                            const polynomial_frobenius_context<Backend> &frobenius_context,
                                            polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if (degree_count == 0 || degree_count > precomputation.block_size()) {
            throw std::invalid_argument("the coarse Kaltofen-Shoup degree count must be within the block");
        }

        polynomial_type interval_product = {value_type::one()};
        polynomial_type difference;
        for (std::size_t offset = 0; offset < degree_count; ++offset) {
            const std::size_t baby_step_index = precomputation.block_size() - 1 - offset;
            subtraction(difference, giant_step, precomputation.baby_step(baby_step_index));
            mulmod(interval_product, interval_product, difference, frobenius_context.divisor_context(),
                   arithmetic_context);
        }
        gcd(output, remaining, interval_product, arithmetic_context);
    }

    /**
     * Split one coarse Kaltofen-Shoup block into exact-degree groups. Let the coefficient field contain Q elements,
     * let B be the modulus used to construct precomputation, let b = precomputation.block_size(), and let
     * h[i] = X^(Q^i) mod B. If giant_step = X^(Q^e) mod B, the first degree in its block is
     *
     *     first_factor_degree = e - b + 1.
     *
     * The polynomial passed as coarse_block contains only irreducible factors whose degrees range from
     * first_factor_degree through first_factor_degree + degree_count - 1. At offset i, subtracting h[b - 1 - i] from
     * giant_step targets degree
     *
     *     e - (b - 1 - i) = first_factor_degree + i.
     *
     * After the lower degrees in this block have been removed, taking the GCD with coarse_block therefore extracts
     * exactly the factors of that degree.
     *
     * The extracted groups are appended in increasing degree order. The callback is invoked immediately after each
     * group is appended and may stop the split early.
     *
     * @return stop_factorization if the callback requests an early stop; continue_factorization otherwise.
     * @throws std::invalid_argument if first_factor_degree is zero, or if degree_count is zero or greater than the
     * block size.
     * @pre coarse_block is monic and square-free, divides the modulus used for precomputation, and contains only
     * factors in the stated degree interval.
     * @pre giant_step and first_factor_degree describe the same block as precomputation.
     */
    template<SupportsDivrem Backend, typename FactorCallback>
        requires DistinctDegreeFactorCallback<FactorCallback, typename Backend::polynomial_type>
    factorization_control kaltofen_shoup_split_coarse_block(
        std::vector<distinct_degree_factor<typename Backend::polynomial_type>> &output,
        typename Backend::polynomial_type coarse_block, const typename Backend::polynomial_type &giant_step,
        std::size_t first_factor_degree, std::size_t degree_count,
        const kaltofen_shoup_frobenius_precomputation<Backend> &precomputation,
        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context, FactorCallback &&factor_callback) {
        using polynomial_type = typename Backend::polynomial_type;

        if (first_factor_degree == 0 || degree_count == 0 || degree_count > precomputation.block_size()) {
            throw std::invalid_argument(
                "the first factor degree and degree count must be positive, and the degree count must not exceed the "
                "block size");
        }

        polynomial_type difference;
        polynomial_type factor;
        polynomial_type quotient;
        for (std::size_t offset = 0; offset < degree_count && coarse_block.size() > 1; ++offset) {
            const std::size_t baby_step_index = precomputation.block_size() - 1 - offset;
            subtraction(difference, giant_step, precomputation.baby_step(baby_step_index));
            gcd(factor, coarse_block, difference, arithmetic_context);
            if (factor.size() == 1) {
                continue;
            }

            output.push_back({std::move(factor), first_factor_degree + offset});
            if (factor_callback(output.back()) == factorization_control::stop_factorization) {
                return factorization_control::stop_factorization;
            }

            factorization_exact_quotient(quotient, coarse_block, output.back().polynomial, arithmetic_context);
            coarse_block = std::move(quotient);
        }

        if (coarse_block.size() > 1) {
            throw std::invalid_argument("the coarse block contains factors outside the stated degree interval");
        }
        return factorization_control::continue_factorization;
    }

    /**
     * Factor a monic square-free polynomial into distinct-degree groups using an explicit Kaltofen-Shoup block size.
     * For each consecutive degree block, the algorithm computes one coarse GCD, removes that block from the
     * unclassified polynomial, and immediately fine-splits it. Interleaving the coarse and fine phases avoids storing
     * a giant-step polynomial for every block while preserving increasing degree order and immediate callback stops.
     *
     * Once the smallest unprocessed factor degree is greater than half the degree of the unclassified polynomial,
     * that polynomial is either constant or irreducible. Otherwise it would contain at least two factors of at least
     * the smallest unprocessed degree, whose combined degree would exceed the polynomial's degree. The irreducible
     * remainder can therefore be emitted directly.
     *
     * @return stop_factorization if the callback requests an early stop; continue_factorization otherwise.
     * @throws std::invalid_argument if block_size is zero.
     * @pre input is monic, square-free, nonzero, and nonconstant.
     */
    template<SupportsDivrem Backend, typename FactorCallback>
        requires DistinctDegreeFactorCallback<FactorCallback, typename Backend::polynomial_type>
    factorization_control kaltofen_shoup_factor_monic_square_free(
        std::vector<distinct_degree_factor<typename Backend::polynomial_type>> &output,
        typename Backend::polynomial_type input, std::size_t block_size,
        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context, FactorCallback &&factor_callback) {
        using polynomial_type = typename Backend::polynomial_type;

        if (block_size == 0) {
            throw std::invalid_argument("the Kaltofen-Shoup block size must be positive");
        }
        // A linear polynomial is already one degree-one factor, so no Frobenius precomputation is needed.
        if (input.size() == 2) {
            output.push_back({std::move(input), 1});
            return factor_callback(output.back());
        }

        polynomial_frobenius_context<Backend> frobenius_context(input, arithmetic_context);
        kaltofen_shoup_frobenius_precomputation<Backend> precomputation(block_size, frobenius_context,
                                                                        arithmetic_context);
        polynomial_type giant_step = precomputation.baby_step(block_size);
        polynomial_type unclassified(std::move(input));

        std::size_t first_factor_degree = 1;
        while (first_factor_degree <= (unclassified.size() - 1) / 2) {
            const std::size_t maximum_degree_to_test = (unclassified.size() - 1) / 2;
            const std::size_t degree_count = std::min(block_size, maximum_degree_to_test - first_factor_degree + 1);

            polynomial_type coarse_block;
            kaltofen_shoup_coarse_block_factor(coarse_block, unclassified, giant_step, degree_count, precomputation,
                                               frobenius_context, arithmetic_context);
            if (coarse_block.size() > 1) {
                polynomial_type quotient;
                factorization_exact_quotient(quotient, unclassified, coarse_block, arithmetic_context);
                unclassified = std::move(quotient);

                if (kaltofen_shoup_split_coarse_block(output, std::move(coarse_block), giant_step, first_factor_degree,
                                                      degree_count, precomputation, arithmetic_context,
                                                      factor_callback) == factorization_control::stop_factorization) {
                    return factorization_control::stop_factorization;
                }
            }

            // Giant-step exponents remain aligned to fixed-size blocks, including when the final tested interval is
            // shorter than block_size.
            first_factor_degree += block_size;
            if (first_factor_degree <= (unclassified.size() - 1) / 2) {
                precomputation.apply_block_frobenius(giant_step, giant_step, frobenius_context, arithmetic_context);
            }
        }

        if (unclassified.size() > 1) {
            const std::size_t irreducible_factor_degree = unclassified.size() - 1;
            output.push_back({std::move(unclassified), irreducible_factor_degree});
            if (factor_callback(output.back()) == factorization_control::stop_factorization) {
                return factorization_control::stop_factorization;
            }
        }
        return factorization_control::continue_factorization;
    }

}    // namespace nil::crypto3::math::detail

namespace nil::crypto3::math {

    /**
     * Split a square-free polynomial into products of irreducible factors of equal degree using the blocked
     * Kaltofen-Shoup distinct-degree algorithm. The input is normalized to monic form while its original leading
     * coefficient is preserved in the result. Zero and constant inputs produce no factors.
     *
     * The block size is selected automatically to balance baby and giant Frobenius steps. Factors are emitted in
     * increasing irreducible-factor degree. After each factor is appended, factor_callback may request an early stop;
     * the stopped result includes that factor and has complete set to false.
     *
     * @throws std::invalid_argument if a nonconstant input is not square-free.
     * @pre input is a nonempty coefficient polynomial.
     */
    template<detail::SupportsDivrem Backend, typename FactorCallback>
        requires detail::DistinctDegreeFactorCallback<FactorCallback, typename Backend::polynomial_type>
    distinct_degree_factorization_result<typename Backend::polynomial_type>
        distinct_degree_factorization_kaltofen_shoup(
            const typename Backend::polynomial_type &input,
            polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
            FactorCallback &&factor_callback) {
        using polynomial_type = typename Backend::polynomial_type;
        using result_type = distinct_degree_factorization_result<polynomial_type>;

        result_type result;
        polynomial_type monic_input;
        if (!detail::prepare_distinct_degree_factorization_input<Backend>(monic_input, result.leading_coefficient,
                                                                          input, arithmetic_context)) {
            return result;
        }

        const std::size_t block_size = detail::kaltofen_shoup_block_size(monic_input.size() - 1);
        if (detail::kaltofen_shoup_factor_monic_square_free(result.factors, std::move(monic_input), block_size,
                                                            arithmetic_context, factor_callback) ==
            factorization_control::stop_factorization) {
            result.complete = false;
        }
        return result;
    }

    /** Compute the complete Kaltofen-Shoup distinct-degree factorization without a staged callback. */
    template<detail::SupportsDivrem Backend>
    distinct_degree_factorization_result<typename Backend::polynomial_type>
        distinct_degree_factorization_kaltofen_shoup(
            const typename Backend::polynomial_type &input,
            polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using factor_type = distinct_degree_factor<typename Backend::polynomial_type>;
        return distinct_degree_factorization_kaltofen_shoup<Backend>(
            input, arithmetic_context,
            [](const factor_type &) { return factorization_control::continue_factorization; });
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_KALTOFEN_SHOUP_DISTINCT_DEGREE_FACTORIZATION_HPP
