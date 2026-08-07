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

#ifndef CRYPTO3_MATH_UNITY_ROOT_HPP
#define CRYPTO3_MATH_UNITY_ROOT_HPP

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <complex>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/number.hpp>

#include <nil/crypto3/algebra/type_traits.hpp>
#include <nil/crypto3/algebra/fields/params.hpp>

#include <nil/crypto3/math/totient.hpp>

namespace nil {
    namespace crypto3 {
        namespace math {
            using namespace boost::multiprecision;

            namespace detail {

                // Return the distinct prime divisors of n in increasing order, discarding repeated factors.
                inline std::vector<std::size_t> distinct_prime_factors(std::size_t n) {
                    std::vector<std::size_t> factors;

                    for (std::size_t divisor = 2; divisor <= n / divisor; ++divisor) {
                        if (n % divisor != 0) {
                            continue;
                        }

                        factors.push_back(divisor);
                        do {
                            n /= divisor;
                        } while (n % divisor == 0);
                    }

                    if (n > 1) {
                        factors.push_back(n);
                    }

                    return factors;
                }

            }    // namespace detail

            /*
             A helper function to RootOfUnity function. This finds a generator for a given
             prime q. Input: BigInteger q which is a prime. Output: A generator of prime q
             */
            template<typename IntegerType>
            static typename std::enable_if<!algebra::is_field_element<IntegerType>::value, IntegerType>::type
                find_generator(const IntegerType &q) {
                std::set<IntegerType> prime_factors;

                IntegerType qm1 = q - IntegerType(1);
                IntegerType qm2 = q - IntegerType(2);
                algebra::prime_factorize<IntegerType>(qm1, prime_factors);
                bool generator_found = false;
                IntegerType gen;
                while (!generator_found) {
                    uint32_t count = 0;

                    // gen = RNG(qm2).ModAdd(IntegerType::ONE, q); //modadd note needed
                    gen = UniformRandomBitGenerator(qm2) + IntegerType(1);

                    for (auto it = prime_factors.begin(); it != prime_factors.end(); ++it) {
                        if (gen.ModExp(qm1 / (*it), q) == IntegerType(1))
                            break;
                        else
                            count++;
                    }
                    if (count == prime_factors.size())
                        generator_found = true;
                }
                return gen;
            }

            /**
             * Finds roots of unity for given input.  Assumes the the input is a power of
             * two.
             *
             * @param m as number which is cyclotomic(in format of int).
             * @param &modulo which is used to find generator.
             *
             * finds roots of unity for given input.  Assumes the the input is a power of two.
             Mostly likely does not give correct results otherwise. input:  m as number
             which is cyclotomic(in format of int), modulo which is used to find generator
             (in format of BigInteger)

             output:  root of unity (in format of BigInteger)
             *
             * @return a root of unity.
             */
            template<typename Backend, boost::multiprecision::expression_template_option ExpressionTemplates>
            boost::multiprecision::number<Backend, ExpressionTemplates>
                unity_root(uint32_t m, const boost::multiprecision::number<Backend, ExpressionTemplates> &modulo) {
                using namespace boost::multiprecision;

                number<Backend, ExpressionTemplates> M(m);

                if ((modulo - number<Backend, ExpressionTemplates>(1) % M) % M != 0) {
                    return {};
                }

                number<backends::modular_adaptor<Backend, backends::modular_params_rt<Backend>>, ExpressionTemplates>
                    gen(find_generator(modulo), modulo), result = boost::multiprecision::pow(gen, (modulo - 1) / M);
                if (result == 1u) {
                    result = unity_root(m, modulo);
                }

                /*
                 * At this point, result contains a primitive root of unity. However,
                 * we want to return the minimum root of unity, to avoid different
                 * crypto contexts having different roots of unity for the same
                 * cyclotomic order and moduli. Therefore, we are going to cycle over
                 * all primitive roots of unity and select the smallest one (minRU).
                 *
                 * To cycle over all primitive roots of unity, we raise the root of
                 * unity in result to all the powers that are co-prime to the
                 * cyclotomic order. In power-of-two cyclotomics, this will be the
                 * set of all odd powers, but here we use a more general routine
                 * to support arbitrary cyclotomics.
                 *
                 */

                number<Backend, ExpressionTemplates> mu = modulo.ComputeMu();
                number<Backend, ExpressionTemplates> x(1);
                x.ModMulEq(result, modulo, mu);
                number<Backend, ExpressionTemplates> minRU(x);
                number<Backend, ExpressionTemplates> curPowIdx(1);
                std::vector<number<Backend, ExpressionTemplates>> coprimes =
                    totient_list<number<Backend, ExpressionTemplates>>(m);
                for (uint32_t i = 0; i < coprimes.size(); i++) {
                    auto nextPowIdx = coprimes[i];
                    number<Backend, ExpressionTemplates> diffPow(nextPowIdx - curPowIdx);
                    for (std::size_t j = 0; j < diffPow; j++) {
                        x.ModMulEq(result, modulo, mu);
                    }
                    if (x < minRU && x != number<Backend, ExpressionTemplates>(1)) {
                        minRU = x;
                    }
                    curPowIdx = nextPowIdx;
                }
                return minRU;
            }

            template<typename FieldType>
            constexpr typename std::enable_if<std::is_same<typename FieldType::value_type, std::complex<double>>::value,
                                              typename FieldType::value_type>::type
                unity_root(const std::size_t n) {
                const double PI = boost::math::constants::pi<double>();

                return typename FieldType::value_type(cos(2 * PI / n), sin(2 * PI / n));
            }

            /**
             * Return a field element of order n. Power-of-two orders use Crypto3's existing radix-2 construction;
             * non-power-of-two orders use the field's multiplicative generator.
             */
            template<typename FieldType>
            constexpr
                typename std::enable_if<!std::is_same<typename FieldType::value_type, std::complex<double>>::value,
                                        typename FieldType::value_type>::type
                unity_root(const std::size_t n) {
                typedef typename FieldType::value_type value_type;

                if (n == 0) {
                    throw std::invalid_argument("unity_root: expected n > 0");
                }

                if (n == 1) {
                    return value_type::one();
                }

                const bool is_power_of_two = (n & (n - 1)) == 0;
                if (!is_power_of_two) {
                    typedef typename FieldType::integral_type integral_type;

                    const integral_type multiplicative_group_order = FieldType::modulus - 1;
                    if (multiplicative_group_order % n != 0) {
                        throw std::invalid_argument("unity_root: n must divide FieldType::modulus - 1");
                    }

                    const integral_type exponent = multiplicative_group_order / n;
                    const value_type generator =
                        value_type(algebra::fields::arithmetic_params<FieldType>::multiplicative_generator);
                    const value_type omega = generator.pow(exponent);

                    if (omega.pow(n) != value_type::one()) {
                        throw std::invalid_argument("unity_root: failed to construct a root of the requested order");
                    }

                    for (const std::size_t prime_factor : detail::distinct_prime_factors(n)) {
                        if (omega.pow(n / prime_factor) == value_type::one()) {
                            throw std::invalid_argument("unity_root: multiplicative generator has insufficient order");
                        }
                    }

                    return omega;
                }

                const std::size_t logn = std::ceil(std::log2(n));

                if (n != (1u << logn)) {
                    throw std::invalid_argument("expected n == (1u << logn)");
                }
                if (logn > algebra::fields::arithmetic_params<FieldType>::s) {
                    throw std::invalid_argument("expected logn <= arithmetic_params<FieldType>::s");
                }

                value_type omega = value_type(algebra::fields::arithmetic_params<FieldType>::root_of_unity);
                for (std::size_t i = algebra::fields::arithmetic_params<FieldType>::s; i > logn; --i) {
                    omega *= omega;
                }

                return omega;
            }
        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MATH_UNITY_ROOT_HPP
