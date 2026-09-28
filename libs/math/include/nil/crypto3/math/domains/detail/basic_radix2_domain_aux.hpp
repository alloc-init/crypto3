//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2024 Dmitrii Tabalin <d.tabalin@nil.foundation>
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

#ifndef CRYPTO3_MATH_BASIC_RADIX2_DOMAIN_AUX_HPP
#define CRYPTO3_MATH_BASIC_RADIX2_DOMAIN_AUX_HPP

#include <algorithm>
#include <bit>
#include <memory>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/algebra/type_traits.hpp>

#include <nil/crypto3/math/algorithms/batch_inverse.hpp>
#include <nil/crypto3/math/algorithms/unity_root.hpp>
#include <nil/crypto3/math/detail/field_utils.hpp>

#include <nil/crypto3/bench/scoped_profiler.hpp>

namespace nil {
    namespace crypto3 {
        namespace math {
            namespace detail {

                /*
                 * Building caches for fft operations
                 */
                template<typename FieldType>
                void create_fft_cache(const std::size_t size,
                                      const typename FieldType::value_type &omega,
                                      std::vector<typename FieldType::value_type> &cache) {
                    typedef typename FieldType::value_type value_type;
                    cache.resize(size, FieldType::value_type::zero());
                    cache[0] = 1;
                    for (std::size_t i = 1; i < size; ++i) {
                        cache[i] = cache[i - 1] * omega;
                    }
                }

                /*
                 * Below we make use of pseudocode from [CLRS 2n Ed, pp. 864].
                 * Also, note that it's the caller's responsibility to multiply by 1/N.
                 */
                template<typename FieldType, typename Range>
                void basic_radix2_fft_cached(Range &a, const std::vector<typename FieldType::value_type> &omega_cache) {
                    typedef typename std::iterator_traits<decltype(std::begin(std::declval<Range>()))>::value_type
                        value_type;
                    BOOST_STATIC_ASSERT(algebra::Field<FieldType>);

                    // It now supports curve elements too, should probably some other assertion about the field type and
                    // value type BOOST_STATIC_ASSERT(std::is_same<typename FieldType::value_type, value_type>::value);

                    const std::size_t n = a.size(), logn = log2(n);
                    if (n != (1u << logn))
                        throw std::invalid_argument("expected n == (1u << logn)");
                    bench::register_fft<FieldType>(logn);

                    // swapping in place (from Storer's book)
                    for (std::size_t k = 0; k < n; ++k) {
                        const std::size_t rk = crypto3::math::detail::bitreverse(k, logn);
                        if (k < rk)
                            std::swap(a[k], a[rk]);
                    }

                    // invariant: m = 2^{s-1}
                    for (std::size_t s = 1, m = 1, inc = n / 2; s <= logn; ++s, m <<= 1, inc >>= 1) {
                        // w_m is 2^s-th root of unity now
                        for (std::size_t k = 0; k < n; k += 2 * m) {
                            for (std::size_t j = 0, idx = 0; j < m; ++j, idx += inc) {
                                value_type t = a[k + j + m];
                                t *= omega_cache[idx];
                                a[k + j + m] = a[k + j];
                                a[k + j + m] -= t;
                                a[k + j] += t;
                            }
                        }
                    }
                }

                /**
                 * Note that it's the caller's responsibility to multiply by 1/N.
                 */
                template<typename FieldType, typename Range>
                void basic_radix2_fft(
                    Range &a, const typename FieldType::value_type &omega,
                    std::shared_ptr<std::vector<typename FieldType::value_type>> omega_cache = nullptr) {

                    if (omega_cache == nullptr) {
                        std::vector<typename FieldType::value_type> omega_powers;
                        create_fft_cache<FieldType>(a.size(), omega, omega_powers);
                        basic_radix2_fft_cached<FieldType>(a, omega_powers);
                    } else {
                        basic_radix2_fft_cached<FieldType>(a, *omega_cache);
                    }
                }

                /**
                 * Compute the Lagrange coefficients and Z(t) from precomputed domain powers.
                 *
                 * @pre omega_powers contains (1, omega, ..., omega^(m-1)) for a primitive m-th root of unity.
                 */
                template<typename FieldType>
                std::vector<typename FieldType::value_type> basic_radix2_evaluate_all_lagrange_polynomials(
                    const std::vector<typename FieldType::value_type> &omega_powers,
                    const typename FieldType::value_type &t,
                    typename FieldType::value_type &vanishing_polynomial_at_t) {
                    typedef typename FieldType::value_type value_type;

                    const std::size_t m = omega_powers.size();
                    if (!std::has_single_bit(m)) {
                        throw std::invalid_argument("expected a nonzero power-of-two domain size");
                    }

                    vanishing_polynomial_at_t = t.pow(m) - value_type::one();
                    if (m == 1) {
                        return std::vector<value_type>(1, value_type::one());
                    }

                    std::vector<value_type> denominators;
                    denominators.reserve(m + 1);
                    for (std::size_t i = 0; i < m; ++i) {
                        const value_type denominator = t - omega_powers[i];
                        if (denominator.is_zero()) {
                            // A domain point selects one basis polynomial; do not try to invert zero.
                            std::vector<value_type> result(m, value_type::zero());
                            result[i] = value_type::one();
                            return result;
                        }
                        denominators.push_back(denominator);
                    }

                    // Include m so normalization by 1/m uses the same single field inversion.
                    denominators.emplace_back(m);
                    auto result = batch_inverse_nonzero(denominators);
                    const value_type scale = vanishing_polynomial_at_t * result.back();
                    result.pop_back();
                    // L_i(t) = Z(t) * omega^i / (m * (t - omega^i)). Reuse the inverse vector for the result.
                    for (std::size_t i = 0; i < m; ++i) {
                        result[i] *= scale * omega_powers[i];
                    }

                    return result;
                }

                /**
                 * Compute the m Lagrange coefficients for S = (1, omega, ..., omega^(m-1)) at t.
                 */
                template<typename FieldType>
                std::vector<typename FieldType::value_type>
                    basic_radix2_evaluate_all_lagrange_polynomials(const std::size_t m,
                                                                   const typename FieldType::value_type &t) {
                    typedef typename FieldType::value_type value_type;

                    if (!std::has_single_bit(m)) {
                        throw std::invalid_argument("expected a nonzero power-of-two domain size");
                    }
                    if (m == 1) {
                        return std::vector<value_type>(1, value_type::one());
                    }

                    std::vector<value_type> omega_powers;
                    create_fft_cache<FieldType>(m, unity_root<FieldType>(m), omega_powers);
                    value_type vanishing_polynomial_at_t;
                    return basic_radix2_evaluate_all_lagrange_polynomials<FieldType>(omega_powers, t,
                                                                                     vanishing_polynomial_at_t);
                }
            }    // namespace detail
        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // ALGEBRA_FFT_BASIC_RADIX2_DOMAIN_AUX_HPP
