//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
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

#ifndef CRYPTO3_MATH_GEOMETRIC_SEQUENCE_DOMAIN_HPP
#define CRYPTO3_MATH_GEOMETRIC_SEQUENCE_DOMAIN_HPP

#include <stdexcept>
#include <vector>

#include <nil/crypto3/math/algorithms/batch_inverse.hpp>
#include <nil/crypto3/math/domains/evaluation_domain.hpp>

#include <nil/crypto3/math/polynomial/basis_change.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>

namespace nil {
    namespace crypto3 {
        namespace math {

            using namespace nil::crypto3::algebra;

            template<typename FieldType, typename ValueType>
            class evaluation_domain;

            /**
             * An exact-size evaluation domain whose points are x_i = r^i for 0 <= i < m.
             *
             * The field-element Lagrange path precomputes the domain points and barycentric weights in O(m) field
             * operations, using one field inversion through batch inversion. Each subsequent evaluation of all m
             * Lagrange basis polynomials at one point takes O(m) field operations and one field inversion. Therefore,
             * after constructing one reusable domain, evaluating the weights at k points takes O(m * k), not
             * quadratic work per point.
             *
             * This complexity guarantee applies to the field-element Lagrange overloads. The legacy transforms and
             * powers-based overload are documented separately below.
             */
            template<typename FieldType, typename ValueType = typename FieldType::value_type>
            class geometric_sequence_domain : public evaluation_domain<FieldType, ValueType> {
                typedef typename FieldType::value_type field_value_type;
                typedef ValueType value_type;

                struct precomputation_type {
                    field_value_type generator;
                    std::vector<field_value_type> geometric_sequence;
                    std::vector<field_value_type> geometric_triangular_sequence;
                    std::vector<field_value_type> barycentric_weights;

                    /*
                     * Here m is the exact domain size: there are m points x_0, ..., x_(m-1), and transforms consume
                     * exactly m coefficients or evaluations.
                     *
                     * Build the reusable field data once, in three linear stages:
                     *   1. Generate r^i and r^(i(i-1)/2), while validating that the domain points are distinct.
                     *   2. Batch-invert a packed denominator vector using one field inversion. Its first entry is r,
                     *      and entry i > 0 is 1 - r^i, so the result supplies both r^-1 and every
                     *      (1 - r^i)^-1 needed below.
                     *   3. Derive the barycentric weights by recurrence.
                     *
                     * The completed object is immutable and constructor-only scratch remains local, so concurrent
                     * Lagrange evaluations share only read-only state.
                     */
                    explicit precomputation_type(const std::size_t m) :
                        generator(fields::arithmetic_params<FieldType>::geometric_generator),
                        geometric_sequence(m, field_value_type::zero()),
                        geometric_triangular_sequence(m, field_value_type::zero()),
                        barycentric_weights(m, field_value_type::zero()) {
                        if (generator.is_zero()) {
                            throw std::invalid_argument("geometric: expected a nonzero geometric generator");
                        }

                        geometric_sequence[0] = field_value_type::one();
                        geometric_triangular_sequence[0] = field_value_type::one();
                        std::vector<field_value_type> denominators(m, field_value_type::zero());
                        denominators[0] = generator;
                        for (std::size_t i = 1; i < m; ++i) {
                            geometric_sequence[i] = geometric_sequence[i - 1] * generator;
                            if (geometric_sequence[i].is_one()) {
                                throw std::invalid_argument(
                                    "geometric: the geometric generator does not define distinct domain points");
                            }
                            geometric_triangular_sequence[i] =
                                geometric_triangular_sequence[i - 1] * geometric_sequence[i - 1];
                            denominators[i] = field_value_type::one() - geometric_sequence[i];
                        }

                        // Inverting the generator and all nonzero (1 - r^i) terms together makes domain construction
                        // use a single field inversion. Every other inverse below follows by recurrence.
                        const std::vector<field_value_type> inverse_denominators = batch_inverse_nonzero(denominators);

                        std::vector<field_value_type> inverse_geometric_sequence(m, field_value_type::zero());
                        inverse_geometric_sequence[0] = field_value_type::one();
                        for (std::size_t i = 1; i < m; ++i) {
                            inverse_geometric_sequence[i] = inverse_geometric_sequence[i - 1] * inverse_denominators[0];
                        }

                        /*
                         * Let Z(X) = product_(j=0)^(m-1) (X - x_j) be the domain's vanishing polynomial and let
                         * Z'(X) be its formal derivative. At a domain point,
                         *
                         *   Z'(x_i) = product_(j != i) (x_i - x_j).
                         *
                         * The barycentric weight w_i is defined as 1 / Z'(x_i). For geometric points these weights
                         * have the linear recurrence
                         *
                         *   w_0 = product_(j=1)^(m-1) (1 - x_j)^-1,
                         *   w_i = -w_(i-1) x_(m-1-i)^-1 (1 - x_(m-i)) / (1 - x_i).
                         *
                         * This is the ratio of the two adjacent derivative products for geometric points and avoids
                         * the quadratic pairwise products used by a generic barycentric-weight construction.
                         */
                        barycentric_weights[0] = field_value_type::one();
                        for (std::size_t i = 1; i < m; ++i) {
                            barycentric_weights[0] *= inverse_denominators[i];
                        }
                        for (std::size_t i = 1; i < m; ++i) {
                            barycentric_weights[i] =
                                -(barycentric_weights[i - 1] * inverse_geometric_sequence[m - 1 - i] *
                                  denominators[m - i] * inverse_denominators[i]);
                        }
                    }
                };

                static std::size_t validate_size(const std::size_t m) {
                    if (m <= 1) {
                        throw std::invalid_argument("geometric: expected m > 1");
                    }
                    return m;
                }

                const precomputation_type precomputation_;

            public:
                typedef FieldType field_type;

                explicit geometric_sequence_domain(const std::size_t m) :
                    evaluation_domain<FieldType, ValueType>(validate_size(m)), precomputation_(m) {
                }

                /*
                 * TODO: These legacy transforms are independent of the field-element Lagrange evaluation below.
                 * They still perform individual inversions and use the generic polynomial-multiplication backend,
                 * which cannot handle the target sizes over low-2-adicity fields. Optimize them separately before
                 * using geometric-domain transforms at those sizes.
                 */
                void fft(std::vector<value_type> &a) override {
                    if (a.size() != this->m) {
                        if (a.size() < this->m) {
                            a.resize(this->m, value_type::zero());
                        } else {
                            throw std::invalid_argument("geometric: expected a.size() == this->m");
                        }
                    }

                    monomial_to_newton_basis_geometric<FieldType>(
                        a, precomputation_.geometric_sequence, precomputation_.geometric_triangular_sequence, this->m);

                    /* Newton to Evaluation */
                    std::vector<field_value_type> T(this->m);
                    T[0] = field_value_type::one();

                    std::vector<value_type> g(this->m);
                    g[0] = a[0];

                    for (std::size_t i = 1; i < this->m; i++) {
                        T[i] = T[i - 1] * (precomputation_.geometric_sequence[i] - field_value_type::one()).inversed();
                        g[i] = precomputation_.geometric_triangular_sequence[i] * a[i];
                    }

                    multiplication(a, g, T);
                    a.resize(this->m);

                    for (std::size_t i = 0; i < a.size(); ++i) {
                        a[i] *= T[i].inversed();
                    }
                }

                void inverse_fft(std::vector<value_type> &a) override {
                    if (a.size() != this->m) {
                        if (a.size() < this->m) {
                            a.resize(this->m, value_type::zero());
                        } else {
                            throw std::invalid_argument("geometric: expected a.size() == this->m");
                        }
                    }

                    /* Interpolation to Newton */
                    std::vector<field_value_type> T(this->m);
                    T[0] = field_value_type::one();

                    std::vector<value_type> W(this->m);
                    W[0] = a[0] * T[0];

                    field_value_type prev_T = T[0];
                    for (std::size_t i = 1; i < this->m; i++) {
                        prev_T *= (precomputation_.geometric_sequence[i] - field_value_type::one()).inversed();

                        W[i] = a[i] * prev_T;
                        T[i] = precomputation_.geometric_triangular_sequence[i] * prev_T;
                        if (i % 2 == 1)
                            T[i] = -T[i];
                    }

                    multiplication(a, W, T);
                    a.resize(this->m);

                    for (std::size_t i = 0; i < a.size(); ++i) {
                        a[i] *= precomputation_.geometric_triangular_sequence[i].inversed();
                    }

                    newton_to_monomial_basis_geometric<FieldType>(
                        a, precomputation_.geometric_sequence, precomputation_.geometric_triangular_sequence, this->m);
                }

                void batch_fft(std::vector<std::vector<value_type>> &a) override {
                    // TODO(martun): implement this.
                    throw std::logic_error {"Not implemented yet"};
                }

                void batch_inverse_fft(std::vector<std::vector<value_type>> &a) override {
                    // TODO(martun): implement this.
                    throw std::logic_error {"Not implemented yet"};
                }

                /**
                 * Evaluate every Lagrange basis polynomial at t and return Z(t) through vanishing_polynomial_at_t.
                 * With
                 *
                 *   Z(X) = product_(j=0)^(m-1) (X - x_j),
                 *   w_i = 1 / Z'(x_i),
                 *
                 * the result away from the domain is
                 *
                 *   L_i(t) = Z(t) w_i / (t - x_i).
                 *
                 * All (t - x_i)^-1 values are obtained with one batch inversion. If t is a domain point, the
                 * corresponding unit vector is returned before attempting inversion.
                 */
                std::vector<field_value_type>
                    evaluate_all_lagrange_polynomials(const field_value_type &t,
                                                      field_value_type &vanishing_polynomial_at_t) override {
                    std::vector<field_value_type> denominators(this->m, field_value_type::zero());
                    vanishing_polynomial_at_t = field_value_type::one();
                    for (std::size_t i = 0; i < this->m; ++i) {
                        denominators[i] = t - precomputation_.geometric_sequence[i];
                        if (denominators[i].is_zero()) {
                            vanishing_polynomial_at_t = field_value_type::zero();
                            std::vector<field_value_type> res(this->m, field_value_type::zero());
                            res[i] = field_value_type::one();
                            return res;
                        }
                        vanishing_polynomial_at_t *= denominators[i];
                    }

                    const std::vector<field_value_type> inverse_denominators = batch_inverse_nonzero(denominators);
                    std::vector<field_value_type> result(this->m, field_value_type::zero());
                    for (std::size_t i = 0; i < this->m; ++i) {
                        result[i] = vanishing_polynomial_at_t * precomputation_.barycentric_weights[i] *
                                    inverse_denominators[i];
                    }

                    return result;
                }

                std::vector<field_value_type> evaluate_all_lagrange_polynomials(const field_value_type &t) override {
                    field_value_type vanishing_polynomial_at_t;
                    return evaluate_all_lagrange_polynomials(t, vanishing_polynomial_at_t);
                }

                std::vector<value_type> evaluate_all_lagrange_polynomials(
                    const typename std::vector<value_type>::const_iterator &t_powers_begin,
                    const typename std::vector<value_type>::const_iterator &t_powers_end) override {
                    if (std::size_t(std::distance(t_powers_begin, t_powers_end)) < this->m) {
                        throw std::invalid_argument(
                            "geometric_sequence_radix2: expected std::distance(t_powers_begin, t_powers_end) >= "
                            "this->m");
                    }

                    /*
                     * TODO: This legacy overload is independent of the optimized field-element path above. It
                     * materializes the Lagrange polynomials and relies on generic polynomial multiplication and
                     * division, so it is unsuitable for large geometric domains. Optimize it separately before
                     * using the powers-based API at those sizes.
                     *
                     * For the m domain points, evaluate each Lagrange basis polynomial using the supplied powers of
                     * t and return the resulting weights.
                     */

                    /* for all i: w[i] = (1 / r) * w[i-1] * (1 - a[i]^m-i+1) / (1 - a[i]^-i) */

                    /**
                     * If t equals one of the geometric progression values,
                     * then output 1 at the right place, and 0 elsewhere.
                     */
                    for (std::size_t i = 0; i < this->m; ++i) {
                        if (precomputation_.geometric_sequence[i] * t_powers_begin[0] == t_powers_begin[1]) {
                            std::vector<value_type> res(this->m, value_type::zero());
                            res[i] = t_powers_begin[0];
                            return res;
                        }
                    }

                    /**
                     * Otherwise, if t does not equal any of the geometric progression values,
                     * then compute each Lagrange coefficient.
                     */
                    std::vector<polynomial<field_value_type>> l(this->m);

                    l[0] =
                        polynomial<field_value_type>({-precomputation_.geometric_sequence[0], field_value_type::one()});

                    std::vector<field_value_type> g(this->m);
                    g[0] = field_value_type::zero();

                    polynomial<field_value_type> l_vanish = l[0];
                    field_value_type g_vanish = field_value_type::one();
                    for (std::size_t i = 1; i < this->m; i++) {
                        l[i] = polynomial<field_value_type>(
                            {-precomputation_.geometric_sequence[i], field_value_type::one()});
                        g[i] = field_value_type::one() - precomputation_.geometric_sequence[i];

                        l_vanish = l_vanish * l[i];
                        g_vanish *= g[i];
                    }

                    field_value_type r = precomputation_.geometric_sequence[this->m - 1].inversed();
                    field_value_type r_i = r;

                    std::vector<field_value_type> g_i(this->m);
                    g_i[0] = g_vanish.inversed();

                    for (std::size_t i = 0; i < this->m; i++) {
                        l[i] = l_vanish / l[i];
                    }

                    std::vector<value_type> result(this->m, value_type::zero());

                    for (std::size_t j = 0; j < l[0].size(); ++j) {
                        result[0] = result[0] + t_powers_begin[j] * l[0][j];
                    }
                    result[0] = result[0] * g_i[0];
                    for (std::size_t i = 1; i < this->m; i++) {
                        g_i[i] = g_i[i - 1] * g[this->m - i] * -g[i].inversed() * precomputation_.geometric_sequence[i];

                        for (std::size_t j = 0; j < l[i].size(); ++j) {
                            result[i] = result[i] + t_powers_begin[j] * l[i][j];
                        }

                        result[i] = result[i] * (r_i * g_i[i]);
                        r_i *= r;
                    }

                    return result;
                }

                // The base interface requires this legacy name. Geometric domains return r, which need not be a root
                // of unity.
                const field_value_type &get_unity_root() override {
                    return precomputation_.generator;
                }

                field_value_type get_domain_element(const std::size_t idx) override {
                    return precomputation_.geometric_sequence[idx];
                }

                field_value_type compute_vanishing_polynomial(const field_value_type &t) override {
                    // Evaluate the domain vanishing polynomial Z(t) = product_(i=0)^(m-1) (t - x_i).
                    field_value_type Z = field_value_type::one();
                    for (std::size_t i = 0; i < this->m; i++) {
                        Z *= (t - precomputation_.geometric_sequence[i]);
                    }
                    return Z;
                }

                polynomial<field_value_type> get_vanishing_polynomial() override {
                    /*
                     * TODO: Generate and cache these coefficients in O(m) with the Gaussian-binomial recurrence.
                     * Repeated polynomial multiplication is quadratic and may request unsupported roots from the
                     * generic multiplication backend over low-2-adicity fields.
                     */
                    polynomial<field_value_type> z({field_value_type::one()});
                    for (std::size_t i = 0; i < this->m; i++) {
                        z = z * polynomial<field_value_type>(
                                    {-precomputation_.geometric_sequence[i], field_value_type::one()});
                    }
                    return z;
                }

                void add_poly_z(const field_value_type &coeff, std::vector<field_value_type> &H) override {
                    if (H.size() != this->m + 1)
                        throw std::invalid_argument("geometric: expected H.size() == this->m+1");

                    // TODO: Reuse the linear-time, cached vanishing-polynomial coefficients once they are available.
                    // Directly multiplying the m linear factors below takes O(m^2) field operations.
                    std::vector<field_value_type> vanishing_polynomial(this->m + 1, field_value_type::zero());
                    vanishing_polynomial[0] = field_value_type::one();
                    for (std::size_t i = 0; i < this->m; ++i) {
                        const field_value_type &point = precomputation_.geometric_sequence[i];
                        for (std::size_t degree = i + 1; degree > 0; --degree) {
                            vanishing_polynomial[degree] =
                                vanishing_polynomial[degree - 1] - point * vanishing_polynomial[degree];
                        }
                        vanishing_polynomial[0] = -(point * vanishing_polynomial[0]);
                    }

                    for (std::size_t i = 0; i < H.size(); ++i) {
                        H[i] += vanishing_polynomial[i] * coeff;
                    }
                }

                void divide_by_z_on_coset(std::vector<field_value_type> &P) override {
                    /*
                     * TODO: Constant scaling is only correct for multiplicative-subgroup domains. At this geometric
                     * domain's coset points g*r^i, entry i must be divided by Z(g*r^i), and a coset intersecting the
                     * domain must be rejected. This legacy method is not used by field-element Lagrange evaluation.
                     */
                    const field_value_type coset =
                        field_value_type(fields::arithmetic_params<FieldType>::multiplicative_generator);
                    const field_value_type Z_inverse_at_coset = compute_vanishing_polynomial(coset).inversed();

                    for (field_value_type &P_i : P) {
                        P_i *= Z_inverse_at_coset;
                    }
                }
            };
        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // ALGEBRA_FFT_GEOMETRIC_SEQUENCE_DOMAIN_HPP
