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

#ifndef CRYPTO3_MATH_EQUAL_DEGREE_FACTORIZATION_HPP
#define CRYPTO3_MATH_EQUAL_DEGREE_FACTORIZATION_HPP

#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/fields/field_order.hpp>

#include <nil/crypto3/math/polynomial/factorization/polynomial_factorization.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_exponentiation.hpp>

namespace nil::crypto3::math {

    namespace detail {

        /**
         * Sample a canonical polynomial whose degree is less than coefficient_count. Each coefficient is obtained
         * directly from the caller-owned generator, so the caller controls both the random source and its seed.
         *
         * Cantor-Zassenhaus splits an equal-degree group G by sampling an element of the quotient ring F[X]/(G).
         * Every quotient-ring element has a unique polynomial representative of degree below degree(G), which is why
         * the algorithm samples this many coefficients. A random scalar is not sufficient because it has the same
         * residue modulo every irreducible factor of G and therefore cannot generally separate those factors.
         * The highest sampled coefficient may be zero: lower-degree and constant representatives are still valid
         * quotient-ring elements. Cantor-Zassenhaus excludes zero and constant representatives before starting a trial
         * because they cannot separate the irreducible factors of G.
         *
         * @throws std::invalid_argument if coefficient_count is zero.
         */
        template<CoefficientPolynomial Polynomial, typename Generator>
            requires std::constructible_from<Polynomial, std::size_t> &&
                     requires(Polynomial &polynomial, Generator &generator) { polynomial[0] = generator(); }
        Polynomial sample_random_polynomial(std::size_t coefficient_count, Generator &generator) {
            if (coefficient_count == 0) {
                throw std::invalid_argument("random polynomial sampling requires a positive coefficient count");
            }

            Polynomial result(coefficient_count);
            for (std::size_t index = 0; index < coefficient_count; ++index) {
                result[index] = generator();
            }
            condense(result);
            return result;
        }

        /** Precompute the exponent shared by all Cantor-Zassenhaus trials for one irreducible-factor degree. */
        template<polynomial_arithmetic::PolynomialBackend Backend>
        class cantor_zassenhaus_context {
        public:
            using polynomial_type = typename Backend::polynomial_type;
            using value_type = typename polynomial_type::value_type;
            using field_type = typename value_type::field_type;

            explicit cantor_zassenhaus_context(std::size_t irreducible_factor_degree) :
                irreducible_factor_degree_(irreducible_factor_degree) {
                if (irreducible_factor_degree_ == 0) {
                    throw std::invalid_argument("Cantor-Zassenhaus splitting requires a positive factor degree");
                }
                if (algebra::fields::field_characteristic<field_type>() == 2) {
                    throw std::invalid_argument(
                        "Cantor-Zassenhaus splitting for characteristic two is not implemented");
                }

                exponent_ = (algebra::fields::extension_field_order<field_type>(irreducible_factor_degree_) - 1) / 2;
            }

            std::size_t irreducible_factor_degree() const {
                return irreducible_factor_degree_;
            }

            const boost::multiprecision::cpp_int &exponent() const {
                return exponent_;
            }

        private:
            std::size_t irreducible_factor_degree_;
            boost::multiprecision::cpp_int exponent_;
        };

        /**
         * Try once to split a product G of distinct irreducible polynomials that all have the same degree,
         * irreducible_factor_degree, using the odd-characteristic Cantor-Zassenhaus algorithm. Zero and constant random
         * candidates are discarded because they cannot separate the factors of G. The trial uses the first remaining
         * polynomial a, whose degree is between one and degree(G) - 1. If gcd(a, G) is already a proper factor, it is
         * returned immediately. Otherwise compute
         *
         *     quadratic_character = a^((Q^irreducible_factor_degree - 1) / 2) mod G,
         *
         * where Q is the coefficient field's order and d = irreducible_factor_degree. For every irreducible degree-d
         * factor Qi, the quotient F[X]/(Qi) is a field with Q^d elements, and its nonzero multiplicative group has
         * order Q^d - 1. Raising a nonzero residue to half that order gives a value whose square is 1, so in odd
         * characteristic it is either 1 or -1. Thus, modulo each irreducible factor of G, quadratic_character is
         * either 1 or -1, and gcd(quadratic_character - 1, G) selects the factors on which it is 1. The trial succeeds
         * when this GCD is neither 1 nor G. A failed trial sets factor to zero and returns false; it does not sample
         * again.
         *
         * cantor_zassenhaus_context supplies the factor degree and precomputed exponent shared by every trial.
         * divisor_context is supplied separately so its polynomial inverse can be reused by subsequent trials against
         * this particular G. Its divisor must be monic and square-free, all its irreducible factors must have the
         * context's irreducible_factor_degree, and it must contain at least two such factors.
         *
         * @throws std::invalid_argument if the divisor degree is inconsistent with the requested factor degree or the
         * divisor is not monic.
         */
        template<SupportsDivrem Backend, typename Generator>
        bool try_cantor_zassenhaus_split(typename Backend::polynomial_type &factor,
                                         const cantor_zassenhaus_context<Backend> &split_context,
                                         const polynomial_divisor_context<Backend> &divisor_context,
                                         polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                         Generator &generator) {
            using polynomial_type = typename Backend::polynomial_type;
            using value_type = typename polynomial_type::value_type;

            const std::size_t irreducible_factor_degree = split_context.irreducible_factor_degree();
            const std::size_t group_degree = divisor_context.degree();
            if (group_degree <= irreducible_factor_degree || group_degree % irreducible_factor_degree != 0) {
                throw std::invalid_argument(
                    "Cantor-Zassenhaus splitting requires at least two factors of the requested degree");
            }
            const polynomial_type &group = divisor_context.divisor();
            if (group.back() != value_type::one()) {
                throw std::invalid_argument("Cantor-Zassenhaus splitting requires a monic polynomial");
            }
            polynomial_type random_polynomial;
            // Zero and constant residues have the same value modulo every irreducible factor, so they cannot split G.
            // Resample them before paying for either GCD or modular exponentiation.
            do {
                random_polynomial = sample_random_polynomial<polynomial_type>(group_degree, generator);
            } while (random_polynomial.size() <= 1);
            gcd(factor, random_polynomial, group, arithmetic_context);
            if (factor.size() > 1 && factor.size() < group.size()) {
                return true;
            }

            polynomial_type quadratic_character;
            powmod(quadratic_character, random_polynomial, split_context.exponent(), divisor_context,
                   arithmetic_context);
            quadratic_character[0] = quadratic_character[0] - value_type::one();
            condense(quadratic_character);
            gcd(factor, quadratic_character, group, arithmetic_context);
            if (factor.size() > 1 && factor.size() < group.size()) {
                return true;
            }

            factor.assign(1, value_type::zero());
            return false;
        }

        /**
         * Repeatedly apply Cantor-Zassenhaus splitting until every returned factor has
         * split_context.irreducible_factor_degree(). An explicit worklist avoids recursion depth proportional to the
         * number of factors. Each composite subgroup gets one divisor context that is reused across all failed trials
         * against that subgroup.
         *
         * A successful trial divides the current subgroup into the returned proper factor and its exact quotient, then
         * places both pieces back on the worklist. A piece whose degree equals the requested irreducible-factor degree
         * is already irreducible under the equal-degree input precondition and is moved to the result.
         *
         * @pre group is monic and square-free, all its irreducible factors have the context's factor degree, and
         * generator supplies independent uniformly distributed coefficient-field elements.
         */
        template<SupportsDivrem Backend, typename Generator, typename FactorCallback>
        factorization_control
            cantor_zassenhaus_split_all(typename Backend::polynomial_type group,
                                        const cantor_zassenhaus_context<Backend> &split_context,
                                        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                        Generator &generator, FactorCallback &&factor_callback) {
            using polynomial_type = typename Backend::polynomial_type;

            std::vector<polynomial_type> pending;
            pending.push_back(std::move(group));

            while (!pending.empty()) {
                polynomial_type current = std::move(pending.back());
                pending.pop_back();

                const std::size_t current_degree = current.size() - 1;
                if (current_degree < split_context.irreducible_factor_degree() ||
                    current_degree % split_context.irreducible_factor_degree() != 0) {
                    throw std::invalid_argument(
                        "Cantor-Zassenhaus splitting requires the factor degree to divide every pending degree");
                }
                if (current_degree == split_context.irreducible_factor_degree()) {
                    if (factor_callback(std::move(current)) == factorization_control::stop_factorization) {
                        return factorization_control::stop_factorization;
                    }
                    continue;
                }

                polynomial_divisor_context<Backend> divisor_context(current, current_degree - 1, arithmetic_context);
                polynomial_type factor;
                // The condition performs one complete random trial. Failure requires no state update other than the
                // generator advancing, so the loop body is intentionally empty and the next condition retries.
                while (!try_cantor_zassenhaus_split<Backend>(factor, split_context, divisor_context, arithmetic_context,
                                                             generator)) {
                }

                polynomial_type quotient;
                factorization_exact_quotient(quotient, current, factor, arithmetic_context);
                pending.push_back(std::move(factor));
                pending.push_back(std::move(quotient));
            }

            return factorization_control::continue_factorization;
        }

        /** Split an entire equal-degree group and collect every irreducible factor. */
        template<SupportsDivrem Backend, typename Generator>
        std::vector<typename Backend::polynomial_type> cantor_zassenhaus_split_all(
            typename Backend::polynomial_type group, const cantor_zassenhaus_context<Backend> &split_context,
            polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context, Generator &generator) {
            using polynomial_type = typename Backend::polynomial_type;

            std::vector<polynomial_type> factors;
            cantor_zassenhaus_split_all<Backend>(std::move(group), split_context, arithmetic_context, generator,
                                                 [&](polynomial_type &&factor) {
                                                     factors.push_back(std::move(factor));
                                                     return factorization_control::continue_factorization;
                                                 });
            return factors;
        }

        /**
         * Factor one group produced by distinct-degree factorization and emit its monic irreducible factors directly to
         * factor_callback. This is the composition boundary used by full factorization: it avoids repeating the
         * normalization and square-free GCD performed by the preceding stages, and it does not build an intermediate
         * factorization result. The full-factorization caller can therefore attach the multiplicity inherited from the
         * square-free stage as each irreducible factor is emitted.
         *
         * The inexpensive structural properties are still checked. Square-freeness and the stated common irreducible
         * degree are trusted because verifying them would repeat the preceding factorization stages.
         *
         * @return stop_factorization if factor_callback requests an early stop; continue_factorization otherwise.
         * @throws std::invalid_argument if the group is constant, nonmonic, has a zero factor degree, or its factor
         * degree does not divide its polynomial degree.
         * @pre group.polynomial is canonical and square-free, and all its irreducible factors have
         * group.irreducible_factor_degree.
         */
        template<SupportsDivrem Backend, typename Generator, typename FactorCallback>
        factorization_control
            factor_distinct_degree_group(distinct_degree_factor<typename Backend::polynomial_type> group,
                                         polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                         Generator &generator, FactorCallback &&factor_callback) {
            using value_type = typename Backend::polynomial_type::value_type;

            if (group.irreducible_factor_degree == 0) {
                throw std::invalid_argument("equal-degree factorization requires a positive factor degree");
            }
            if (group.polynomial.size() <= 1) {
                throw std::invalid_argument("equal-degree factorization requires a nonconstant degree group");
            }
            if (group.polynomial.back() != value_type::one()) {
                throw std::invalid_argument("equal-degree factorization requires a monic degree group");
            }
            if ((group.polynomial.size() - 1) % group.irreducible_factor_degree != 0) {
                throw std::invalid_argument(
                    "equal-degree factorization requires the factor degree to divide the polynomial degree");
            }

            cantor_zassenhaus_context<Backend> split_context(group.irreducible_factor_degree);
            return cantor_zassenhaus_split_all<Backend>(std::move(group.polynomial), split_context, arithmetic_context,
                                                        generator, std::forward<FactorCallback>(factor_callback));
        }

        /**
         * Prepare one group produced by distinct-degree factorization for equal-degree splitting. The input is
         * normalized to monic form and its original leading coefficient is returned separately. Constant inputs
         * return false because they contain no factors.
         *
         * Equal-degree factorization requires a square-free input whose irreducible factors all have
         * irreducible_factor_degree. This function checks square-freeness and the necessary divisibility of the total
         * degree. The caller must supply a group produced by distinct-degree factorization. Verifying the degree of
         * every irreducible factor here would require repeating that factorization, so this function does not perform
         * that check.
         *
         * For example, suppose distinct-degree factorization produces the group G = Q1 * Q2 * Q3, where Q1, Q2,
         * and Q3 are distinct irreducible quadratic polynomials. Equal-degree factorization of G with
         * irreducible_factor_degree = 2 separates that group into Q1, Q2, and Q3. Splitting stops when a resulting
         * piece has degree two: under the equal-degree precondition, that piece is one of Q1, Q2, or Q3 and is already
         * irreducible.
         *
         * @throws std::invalid_argument if irreducible_factor_degree is zero, the input is not square-free, or its
         * total degree is not divisible by irreducible_factor_degree.
         * @pre input is a nonempty coefficient polynomial.
         */
        template<SupportsDivrem Backend>
        bool prepare_equal_degree_factorization_input(
            typename Backend::polynomial_type &monic_input,
            typename Backend::polynomial_type::value_type &leading_coefficient,
            const typename Backend::polynomial_type &input, std::size_t irreducible_factor_degree,
            polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            if (irreducible_factor_degree == 0) {
                throw std::invalid_argument("equal-degree factorization requires a positive factor degree");
            }
            if (!prepare_square_free_factorization_input<Backend>(monic_input, leading_coefficient, input,
                                                                  arithmetic_context)) {
                return false;
            }
            if ((monic_input.size() - 1) % irreducible_factor_degree != 0) {
                throw std::invalid_argument(
                    "equal-degree factorization requires the factor degree to divide the polynomial degree");
            }
            return true;
        }

    }    // namespace detail

    /**
     * Split one square-free distinct-degree group into its monic irreducible factors using the odd-characteristic
     * Cantor-Zassenhaus algorithm. Every irreducible input factor must have irreducible_factor_degree; callers normally
     * obtain this pair from distinct-degree factorization. The necessary total-degree divisibility is checked, but the
     * function does not repeat distinct-degree factorization to verify this precondition.
     *
     * The input is normalized to monic form and its original leading coefficient is preserved in the result. Zero and
     * constant inputs return that coefficient and no factors. Complete results contain every irreducible factor with
     * multiplicity one and reconstruct the input using polynomial_factorization_result's usual convention.
     *
     * generator must return independent uniformly distributed coefficient-field elements. After each factor is added
     * to the result, factor_callback may request an early stop. A stopped result contains that factor, has complete set
     * to false, and does not continue splitting pending subgroups.
     *
     * @throws std::invalid_argument if the factor degree is zero, a nonconstant input is not square-free, its degree is
     * not divisible by the factor degree, or the coefficient field has characteristic two.
     * @pre input is a nonempty coefficient polynomial whose irreducible factors all have
     * irreducible_factor_degree.
     */
    template<detail::SupportsDivrem Backend, typename Generator, typename FactorCallback>
        requires detail::PolynomialFactorCallback<FactorCallback, typename Backend::polynomial_type>
    polynomial_factorization_result<typename Backend::polynomial_type>
        equal_degree_factorization(const typename Backend::polynomial_type &input,
                                   std::size_t irreducible_factor_degree,
                                   polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context,
                                   Generator &generator, FactorCallback &&factor_callback) {
        using polynomial_type = typename Backend::polynomial_type;
        using result_type = polynomial_factorization_result<polynomial_type>;

        result_type result;
        polynomial_type monic_input;
        if (!detail::prepare_equal_degree_factorization_input<Backend>(monic_input, result.leading_coefficient, input,
                                                                       irreducible_factor_degree, arithmetic_context)) {
            return result;
        }

        detail::cantor_zassenhaus_context<Backend> split_context(irreducible_factor_degree);
        const factorization_control control = detail::cantor_zassenhaus_split_all<Backend>(
            std::move(monic_input), split_context, arithmetic_context, generator, [&](polynomial_type &&factor) {
                result.factors.push_back({std::move(factor), 1});
                return factor_callback(result.factors.back());
            });
        if (control == factorization_control::stop_factorization) {
            result.complete = false;
        }
        return result;
    }

    /** Compute the complete equal-degree factorization without a staged callback. */
    template<detail::SupportsDivrem Backend, typename Generator>
    polynomial_factorization_result<typename Backend::polynomial_type> equal_degree_factorization(
        const typename Backend::polynomial_type &input, std::size_t irreducible_factor_degree,
        polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context, Generator &generator) {
        using factor_type = polynomial_factor<typename Backend::polynomial_type>;
        return equal_degree_factorization<Backend>(
            input, irreducible_factor_degree, arithmetic_context, generator,
            [](const factor_type &) { return factorization_control::continue_factorization; });
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_EQUAL_DEGREE_FACTORIZATION_HPP
