//---------------------------------------------------------------------------//
// Copyright (c) 2024 Mikhail Komarov <nemo@allocin.it>
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

#ifndef CRYPTO3_ALGEBRA_DLOG_H
#define CRYPTO3_ALGEBRA_DLOG_H

#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/miller_rabin.hpp>

#include <nil/crypto3/algebra/fields/params.hpp>

#include <nil/crypto3/multiprecision/inverse.hpp>
#include <nil/crypto3/functional/internal/errors.hpp>

/**
 * \file
 * \brief Algorithms for computing discrete logarithms.
 *
 * FE schemes instantiated from the Discrete Diffie-Hellman assumption (DDH)
 * all rely on efficient algorithms for calculating discrete logarithms.
 */

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace detail {
                template<typename FieldType>
                void iterate(typename FieldType::value_type &x, typename FieldType::value_type &a,
                             typename FieldType::value_type &b, const typename FieldType::value_type &h) {
                    if (x % 3 == 0) {
                        x = (x * x) % FieldType::modulus;
                        a = (a * 2) % FieldType::modulus;
                        b = (b * 2) % FieldType::modulus;
                    } else if (x % 3 == 1) {
                        x = (x * fields::arithmetic_params<FieldType>::arithmetic_generator) % FieldType::modulus;
                        a = (a + 1) % FieldType::modulus;
                    } else {
                        x = (x * h) % FieldType::modulus;
                        b = (b + 1) % FieldType::modulus;
                    }
                }
            }    // namespace detail

            /**
             * @brief Baby-step giant-step method for computing the discrete logarithm in
             * the Zp group.
             *
             * It searches for a solution <= bound. If bound argument is nil,
             * the bound is automatically set to p-1.
             * The function returns x, where h = g^x mod p. If the solution was not found
             * within the provided bound, it returns an error.
             *
             * @param res Discrete logarithm (the result value placeholder)
             * @param h Element
             * @param g Generator
             * @param p Modulus
             * @param order Order
             * @param bound Bound for solution
             * @return Error code
             */
            template<typename FieldType>
            typename FieldType::value_type baby_giant_dlog(const typename FieldType::value_type &h,
                                                           const typename FieldType::value_type &g,
                                                           const typename FieldType::value_type &bound = 0) {
                using value_type = typename FieldType::value_type;

                // Access modulus and generator from FieldType
                const value_type &p = FieldType::modulus;    // Modulus

                // Check if p is prime
                if (!boost::multiprecision::miller_rabin_test(p, 25)) {
                    return -1;
                }

                // Calculate m = ceil(sqrt(bound or p - 1))
                value_type m =
                    bound != 0 ? boost::multiprecision::sqrt(bound) + 1 : boost::multiprecision::sqrt(p - 1) + 1;

                value_type x = 1;

                // Hash table for the baby-step
                std::unordered_map<value_type, value_type> T;

                // Baby-step: Compute T[g^j mod p] = j for j in [0, m)
                for (value_type j = 0; j < m; ++j) {
                    T[x] = j;
                    x = (x * g) % p;    // x = g^j mod p
                }

                // Compute z = g^(-m) mod p
                value_type z = boost::multiprecision::inverse_mod(g, p);
                z = boost::multiprecision::powm(z, m, p);

                // Giant-step: Check for h * g^(-i*m) mod p in T
                x = h;    // Start with h
                for (value_type i = 0; i < m; ++i) {
                    auto it = T.find(x);
                    if (it != T.end()) {
                        // Found: Compute result
                        return (i * m + it->second) % (p - 1);    // Implicit modular arithmetic
                    }
                    x = (x * z) % p;    // x = h * g^(-i*m) mod p
                }

                return -2;    // Logarithm not found
            }

            template<typename FieldType>
            typename FieldType::value_type inline baby_giant_dlog(const typename FieldType::value_type &h,
                                                                  const typename FieldType::value_type &bound = 0) {
                return baby_giant_dlog<FieldType>(h, fields::arithmetic_params<FieldType>::arithmetic_generator, bound);
            }

            /**
             * @brief Baby-step giant-step method for computing the discrete logarithm in
             * the Zp group finding also negative solutions.
             *
             * It searches for a solution (-bound, bound). If bound argument is nil,
             * the bound is automatically set to p-1 and it works identically than
             * function baby_step_giant_step.
             * The function returns x, where h = g^x mod p. If the solution was not found
             * within the provided bound, it returns an error.
             *
             * @param res Discrete logarithm (the result value placeholder)
             * @param h Element
             * @param g Generator
             * @param p Modulus
             * @param _order Order
             * @param bound Bound for solution
             * @return Error code
             */
            template<typename FieldType>
            typename FieldType::value_type baby_giant_dlog_with_neg(const typename FieldType::value_type &h,
                                                                    const typename FieldType::value_type bound = 0) {
                using value_type = typename FieldType::value_type;

                // Attempt to compute the discrete logarithm with the original generator
                typename FieldType::value_type res = baby_giant_dlog<FieldType>(h, bound);

                if (res != -2 && res != -1) {
                    return res;
                }

                // Compute the inverse of the generator
                const value_type &p = FieldType::modulus;
                const value_type &g = algebra::fields::arithmetic_params<FieldType>::arithmetic_generator;
                value_type g_inv = boost::multiprecision::inverse_mod(g, p);

                // Attempt to compute the discrete logarithm with the inverse generator
                res = baby_giant_dlog<FieldType>(res, h, bound);
                if (res != -2 && res != -1) {
                    res = (p - res) % p;    // Negate the result (equivalent to -res mod p)
                }

                return res;
            }

            /**
             * @brief Baby-step giant-step method for computing the discrete logarithm in
             * the pairing group FP12_BN254 finding also negative solutions.
             *
             * It searches for a solution (-bound, bound). The function returns x, where
             * h = g^x in the group. If the solution was not found within the provided
             * bound, it returns an error.
             *
             * @param res Discrete logarithm (the result value placeholder)
             * @param h Element
             * @param g Generator
             * @param bound Bound for solution
             * @return Error code
             */
            template<typename FieldType>
            typename FieldType::value_type baby_giant_dlog_with_neg(const typename FieldType::value_type &h,
                                                                    const boost::multiprecision::cpp_int &bound) {
                using value_type = typename FieldType::value_type;
                using HashTable = std::unordered_map<value_type, boost::multiprecision::cpp_int>;

                value_type one = FieldType::value_type::one();    // Multiplicative identity

                // Check if h is the identity element
                if (h == one) {
                    return FieldType::value_type::zero();    // Set result to 0
                }

                typename value_type::integral_type m = boost::multiprecision::sqrt(bound) + 1;

                // Baby-step hash table
                HashTable T;

                value_type x = one;

                // Baby-step: Compute T[x] = i for i in [0, m)
                for (typename value_type::integral_type i = 0; i <= m; ++i) {
                    T[x] = i;
                    x = x * algebra::fields::arithmetic_params<FieldType>::arithmetic_generator;    // x = x * g
                }

                // Precompute z = g^(-m)
                value_type z = algebra::fields::arithmetic_params<FieldType>::arithmetic_generator.inverse();
                z = z ^ m;

                // Giant-step: Simultaneously check positive and negative values
                value_type x_neg = h.inverse();    // x_neg = h^(-1)

                x = h;    // Start with h
                for (typename value_type::integral_type i = 0; i <= m; ++i) {
                    // Check x
                    auto it = T.find(x);
                    if (it != T.end()) {
                        return i * m + it->second;    // Compute result
                    }

                    // Check x_neg
                    it = T.find(x_neg);
                    if (it != T.end()) {
                        return -(i * m + it->second);    // Compute negative result
                    }

                    // Update x and x_neg
                    x = x * z;
                    x_neg = x_neg * z;
                }

                return -1;
            }

            /**
             * @brief  Pollard's rho algorithm - simple, non-parallel version.
             *
             * @param res Discrete logarithm (the result value placeholder)
             * @param h Element
             * @param g Generator
             * @param p Modulus
             * @param n Order
             * @return Error code
             */
            template<typename FieldType>
            typename FieldType::value_type pollard_rho_dlog(const typename FieldType::value_type &h) {
                using value_type = typename FieldType::value_type;

                value_type x1 = value_type::one();     // x1 = 1
                value_type a1 = value_type::zero();    // a1 = 0
                value_type b1 = value_type::zero();    // b1 = 0

                value_type x2 = value_type::one();     // x2 = 1
                value_type a2 = value_type::zero();    // a2 = 0
                value_type b2 = value_type::zero();    // b2 = 0

                uint64_t iterations = uint64_t(1) << 32;

                for (uint64_t i = 0; i < iterations; ++i) {
                    // Update x1, a1, b1
                    detail::iterate<FieldType>(x1, a1, b1, h);

                    // Update x2, a2, b2 twice
                    detail::iterate<FieldType>(x2, a2, b2, h);
                    detail::iterate<FieldType>(x2, a2, b2, h);

                    // Check for collision
                    if (x1 == x2) {
                        value_type r = (b2 - b1) % FieldType::modulus;    // b2 - b1 mod n
                        value_type t = (a1 - a2) % FieldType::modulus;    // a1 - a2 mod n

                        if (r == FieldType::value_type(0)) {
                            break;    // Failure: r = 0
                        }

                        // Compute gcd(r, modulus)
                        typename value_type::integral_type gcd_r_modulus =
                            boost::multiprecision::gcd(r, FieldType::modulus);

                        if (gcd_r_modulus == 1) {
                            // Invert r and compute the result
                            return (r.inverse() * t) % FieldType::modulus;
                        }
                        // Handle non-coprime case
                        value_type d = gcd_r_modulus;
                        value_type r_div_d = r / d;
                        value_type t_div_d = t / d;
                        value_type n_div_d = FieldType::modulus / d;

                        value_type inv_r = r_div_d.inverse();
                        value_type q = (inv_r * t_div_d) % n_div_d;

                        for (typename value_type::integral_type j = 0; j < d; ++j) {
                            if ((fields::arithmetic_params<FieldType>::arithmetic_generator ^ q) == h) {
                                return q;
                            }
                            q = (q + n_div_d) % FieldType::modulus;
                        }
                    }
                }

                return -1;
            }
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil

#endif
