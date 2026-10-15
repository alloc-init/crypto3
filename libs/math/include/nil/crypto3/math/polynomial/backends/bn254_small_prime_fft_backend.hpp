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

#ifndef CRYPTO3_MATH_BN254_SMALL_PRIME_FFT_BACKEND_HPP
#define CRYPTO3_MATH_BN254_SMALL_PRIME_FFT_BACKEND_HPP

#include <algorithm>
#include <bit>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/math/detail/bn254_small_prime_fft.hpp>
#include <nil/crypto3/math/polynomial/backends/polynomial_backend.hpp>
#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>

namespace nil::crypto3::math::polynomial_arithmetic {

    /**
     * Exact BN254 Fp12 polynomial products through small-prime transforms and signed CRT.
     * Configure a power-of-two maximum transform length. All smaller power-of-two plans are
     * cached at construction; each operation uses the smallest length containing its full linear
     * product. multiply_low first discards input coefficients that cannot affect the output.
     *
     * One catalogue prefix and coefficient bound cover every supported product, allowing all
     * lengths to share the converter, CRT cofactors and pointwise arithmetic. If a+b-1 <= Nmax,
     * then min(a,b) <= max(1,Nmax/2), so C=102*max(1,Nmax/2)*(p-1)^2 is sufficient. Setup performs
     * the capacity and memory-size checks once; repeated operations require no CRT/FFT setup.
     *
     * Spectrum buffers retain capacity across operations. Squaring converts/transforms one operand
     * and uses the dedicated tower square. Results are assembled separately before swapping into
     * output, permitting output to alias either input and preserving it on failure.
     * Reuse sequentially; concurrent operations require separate backend instances. No constant-time
     * guarantee is provided, and the thread's rounding mode must stay unchanged during an operation.
     */
    class bn254_small_prime_fft_backend {
        using params_type = math::detail::small_prime_fft_params;

    public:
        using field_type = math::detail::bn254_small_prime_fft_converter::field_type;
        using value_type = field_type::value_type;
        using polynomial_type = math::polynomial<value_type>;

        /**
         * max_spectrum_bytes limits one spectrum at the configured maximum length. A second
         * spectrum, cached twiddle tables and coefficient output storage are additional allocations.
         * Invalid/unsupported lengths or an insufficient byte limit throw before allocating plans.
         */
        explicit bn254_small_prime_fft_backend(
            std::size_t max_transform_size,
            std::size_t max_spectrum_bytes = std::numeric_limits<std::size_t>::max()) :
            bn254_small_prime_fft_backend(checked_params(max_transform_size, max_spectrum_bytes)) {
        }

        std::size_t max_transform_size() const {
            return max_transform_size_;
        }

        std::size_t prime_count() const {
            return converter_.prime_count();
        }

        void multiply(polynomial_type &output, const polynomial_type &left, const polynomial_type &right) {
            multiply_low(output, left, right, std::numeric_limits<std::size_t>::max());
        }

        void square(polynomial_type &output, const polynomial_type &input) {
            multiply_prefixes<true>(output, input, input, std::numeric_limits<std::size_t>::max());
        }

        void multiply_low(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                          std::size_t coefficient_count) {
            if (&left == &right) {
                multiply_prefixes<true>(output, left, right, coefficient_count);
            } else {
                multiply_prefixes<false>(output, left, right, coefficient_count);
            }
        }

    private:
        static params_type checked_params(std::size_t max_size, std::size_t max_spectrum_bytes) {
            if (!std::has_single_bit(max_size)) {
                throw std::invalid_argument("bn254_small_prime_fft_backend: expected a power-of-two maximum size");
            }
            // These balanced lengths have product length max_size and attain the maximum term count.
            auto params = math::detail::make_bn254_small_prime_fft_params(
                std::max(std::size_t(1), max_size / 2), max_size / 2 + 1, std::numeric_limits<std::size_t>::max(),
                max_spectrum_bytes);
            if (!params) {
                throw std::invalid_argument(
                    "bn254_small_prime_fft_backend: insufficient roots, CRT capacity or byte limit");
            }
            return std::move(*params);
        }

        explicit bn254_small_prime_fft_backend(const params_type &params) :
            max_transform_size_(params.transform_size), converter_(params.crt.prime_count),
            reconstructor_(params.crt, params.coefficient_bound) {
            const unsigned depth = std::countr_zero(max_transform_size_);
            plans_.reserve((depth + 1) * prime_count());
            pointwise_.reserve(prime_count());
            for (std::size_t prime = 0; prime < prime_count(); ++prime) {
                pointwise_.emplace_back(math::detail::small_prime_fft_primes[prime].modulus);
            }
            for (unsigned level = 0; level <= depth; ++level) {
                for (std::size_t prime = 0; prime < prime_count(); ++prime) {
                    plans_.emplace_back(prime, std::size_t(1) << level);
                }
            }
        }

        template<bool Square>
        void multiply_prefixes(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                               std::size_t coefficient_count) {
            if (plans_.empty()) {
                throw std::invalid_argument("bn254_small_prime_fft_backend: moved-from backend");
            }
            const std::size_t left_size = std::min(left.size(), coefficient_count);
            const std::size_t right_size = std::min(right.size(), coefficient_count);
            if (math::is_zero(left.begin(), left.begin() + left_size) ||
                math::is_zero(right.begin(), right.begin() + right_size)) {
                output.assign(1, value_type::zero());
                return;
            }
            // Nonzero prefixes are nonempty. Check capacity without overflowing left_size+right_size-1.
            if (left_size > max_transform_size_ || right_size > max_transform_size_ - (left_size - 1)) {
                throw std::invalid_argument("bn254_small_prime_fft_backend: product exceeds configured transform size");
            }
            const std::size_t product_size = left_size + right_size - 1;
            const std::size_t size = std::bit_ceil(product_size);
            const std::size_t result_size = std::min(product_size, coefficient_count);
            const std::size_t channels = field_type::arity * prime_count();
            // Constructor checks guarantee that even the largest spectrum's byte count fits size_t.
            left_spectrum_.resize(channels * size);
            converter_.convert(left_spectrum_, std::span<const value_type>(left.data(), left_size), size);
            if constexpr (!Square) {
                right_spectrum_.resize(channels * size);
                converter_.convert(right_spectrum_, std::span<const value_type>(right.data(), right_size), size);
            }
            const std::size_t plan_offset = std::countr_zero(size) * prime_count();
            for (std::size_t prime = 0; prime < prime_count(); ++prime) {
                const auto &plan = plans_[plan_offset + prime];
                auto a = std::span<double>(left_spectrum_)
                             .subspan(prime * field_type::arity * size, field_type::arity * size);
                std::span<double> b;
                if constexpr (!Square) {
                    b = std::span<double>(right_spectrum_)
                            .subspan(prime * field_type::arity * size, field_type::arity * size);
                }
                for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                    plan.fft(a.subspan(coordinate * size, size));
                    if constexpr (!Square) {
                        plan.fft(b.subspan(coordinate * size, size));
                    }
                }
                if constexpr (Square) {
                    pointwise_[prime].square(a, size, plan.size_inverse());
                } else {
                    pointwise_[prime].multiply(a, b, size, plan.size_inverse());
                }
                for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                    plan.inverse_fft(a.subspan(coordinate * size, size), false);
                }
            }
            result_.resize(result_size);
            reconstructor_.reconstruct(std::span<value_type>(result_.data(), result_size), left_spectrum_, size);
            math::condense(result_);
            output.get_storage().swap(result_.get_storage());
        }

        std::size_t max_transform_size_;
        math::detail::bn254_small_prime_fft_converter converter_;
        math::detail::bn254_small_prime_fft_reconstructor reconstructor_;
        std::vector<math::detail::small_prime_fft_plan> plans_;
        std::vector<math::detail::bn254_small_prime_fft_pointwise> pointwise_;
        std::vector<double> left_spectrum_;
        std::vector<double> right_spectrum_;
        polynomial_type result_;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_BN254_SMALL_PRIME_FFT_BACKEND_HPP
