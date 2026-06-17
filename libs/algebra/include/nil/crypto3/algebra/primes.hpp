//---------------------------------------------------------------------------//
// Copyright (c) 2023 Mikhail Komarov <nemo@nil.foundation>
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

#ifndef CRYPTO3_ALGEBRA_PRIMES_HPP
#define CRYPTO3_ALGEBRA_PRIMES_HPP

#include <cassert>

#include <set>

#include <boost/multiprecision/miller_rabin.hpp>

#include <nil/crypto3/algebra/random_element.hpp>

#include <nil/crypto3/multiprecision/modular/modular_adaptor.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            using namespace boost::multiprecision;

            /*
             The Pollard Rho factorization of a number n.
             Input: n the number to be factorized.
             Output: a factor of n.
             */
            template<typename Backend, expression_template_option ExpressionTemplates>
            number<Backend, ExpressionTemplates>
                pollard_rho_factorization(const number<Backend, ExpressionTemplates> &n) {
                if (!(n % 2)) {
                    return 2;
                }

                boost::random::independent_bits_engine<std::mt19937, 256, number<Backend, ExpressionTemplates>> rng;
                number<backends::modular_adaptor<Backend, backends::modular_params_rt<Backend>>, ExpressionTemplates>
                    divisor, c(rng(), n), x(rng(), n), nn(n, n), xx = x;
                do {
                    x = x * x + c;
                    xx = xx * xx + c;
                    xx = xx * xx + c;
                    divisor = gcd((x > xx) ? x - xx : xx - x, nn);
                } while (static_cast<int>(divisor) == 1);
                return static_cast<number<Backend, ExpressionTemplates>>(divisor);
            }

            /*
             Recursively factorizes and find the distinct primefactors of a number
             Input: n is the number to be prime factorized,
             prime_factors is a set of prime factors of n.
             */
            template<typename IntegerType, std::size_t Iterations = 100>
            void prime_factorize(IntegerType n, std::set<IntegerType> &prime_factors) {
                if (n == 0 || n == 1)
                    return;
                if (miller_rabin_test(n, Iterations)) {
                    prime_factors.insert(n);
                    return;
                }

                IntegerType divisor(pollard_rho_factorization(n));
                IntegerType n_div = n / divisor;
                prime_factorize(divisor, prime_factors);
                prime_factorize(n_div, prime_factors);
            }

            template<typename IntegerType, std::size_t Iterations = 100>
            IntegerType first_prime(uint64_t bits, uint64_t m) {
                IntegerType mi(m), q_new(IntegerType(1) << bits), r(q_new % mi), q_new2(q_new + 1);
                if (r > IntegerType(0))
                    q_new2 += (mi - r);
                BOOST_ASSERT_MSG(q_new2 >= q_new, "FirstPrime parameters overflow this integer implementation");
                while (!miller_rabin_test((q_new = q_new2), Iterations)) {
                    q_new2 = q_new + mi;
                    BOOST_ASSERT_MSG(q_new2 >= q_new, "FirstPrime overflow growing candidate");
                }
                return q_new;
            }

            template<typename IntegerType, std::size_t Iterations = 100>
            IntegerType next_prime(const IntegerType &q, uint64_t m) {
                IntegerType M(m), q_new(q + M);
                while (!miller_rabin_test(q_new, Iterations)) {
                    BOOST_VERIFY_MSG((q_new += M) >= q, "NextPrime overflow growing candidate");
                }
                return q_new;
            }

            template<typename IntegerType, std::size_t Iterations = 100>
            IntegerType previous_prime(const IntegerType &q, uint64_t m) {
                IntegerType M(m), q_new(q - M);
                while (!miller_rabin_test(q_new, Iterations)) {
                    BOOST_VERIFY_MSG((q_new -= M) <= q, "Moduli size is not sufficient! Must be increased.");
                }
                return q_new;
            }
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_PRIMES_HPP
