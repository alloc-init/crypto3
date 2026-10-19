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
#include <array>
#include <bit>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/math/detail/bn254_small_prime_fft.hpp>
#include <nil/crypto3/math/polynomial/backends/polynomial_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>

namespace nil::crypto3::math::polynomial_arithmetic {

    /**
     * Exact BN254 Fp12 polynomial products through small-prime transforms and signed CRT.
     * Configure a power-of-two maximum transform length. All smaller power-of-two plans are
     * cached at construction; ordinary products use the smallest length containing their full linear
     * product. Scalar and linear operands use direct coefficient arithmetic. multiply_low first
     * discards input coefficients that cannot affect the output. Prepared
     * cyclic remainder reconstruction can use a shorter transform once the exact quotient is known.
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
         * Owned, immutable spectrum of a fixed operand's low prefix. The preparation has its own
         * allocation and survives changes to the source polynomial or the preparing backend.
         * Nonzero preparations identify the transform length and catalogue prefix; compatible
         * backend instances can reuse them without retaining references to plans or workspace.
         */
        class prepared_low_product_type {
        public:
            std::size_t coefficient_count() const noexcept {
                return coefficient_count_;
            }

            // Zero denotes a zero fixed prefix or zero requested precision, requiring no spectrum.
            std::size_t transform_size() const noexcept {
                return transform_size_;
            }

            std::size_t storage_bytes() const noexcept {
                return spectrum_.capacity() * sizeof(double);
            }

        private:
            friend class bn254_small_prime_fft_backend;

            prepared_low_product_type(std::size_t coefficient_count, std::size_t prefix_size,
                                      std::size_t transform_size = 0, std::size_t prime_count = 0,
                                      std::vector<double> spectrum = {}) :
                coefficient_count_(coefficient_count), prefix_size_(prefix_size), transform_size_(transform_size),
                prime_count_(prime_count), spectrum_(std::move(spectrum)) {
            }

            std::size_t coefficient_count_;
            std::size_t prefix_size_;
            std::size_t transform_size_;
            std::size_t prime_count_;
            std::vector<double> spectrum_;
        };

        // Owns the transformed complete divisor modulo X^N - 1, including its leading coefficient.
        // Its distinct type keeps cyclic reconstruction separate from multiplication modulo X^k.
        class prepared_cyclic_remainder_type {
        public:
            std::size_t divisor_degree() const noexcept {
                return divisor_degree_;
            }

            std::size_t transform_size() const noexcept {
                return transform_size_;
            }

            std::size_t storage_bytes() const noexcept {
                return spectrum_.capacity() * sizeof(double);
            }

        private:
            friend class bn254_small_prime_fft_backend;

            prepared_cyclic_remainder_type(std::size_t divisor_degree, std::size_t transform_size,
                                           std::size_t prime_count, std::vector<double> spectrum) :
                divisor_degree_(divisor_degree), transform_size_(transform_size), prime_count_(prime_count),
                spectrum_(std::move(spectrum)) {
            }

            std::size_t divisor_degree_;
            std::size_t transform_size_;
            std::size_t prime_count_;
            std::vector<double> spectrum_;
        };

        /**
         * max_spectrum_bytes limits one spectrum at the configured maximum length. A second
         * spectrum, cached twiddle tables and coefficient output storage are additional allocations.
         * Invalid/unsupported lengths or an insufficient byte limit throw before allocating plans.
         * max_matrix_spectrum_bytes limits the additional shared input spectra used by matrix
         * products. Zero disables that optional path; exceeding its budget declines the operation.
         */
        explicit bn254_small_prime_fft_backend(std::size_t max_transform_size,
                                               std::size_t max_spectrum_bytes = std::numeric_limits<std::size_t>::max(),
                                               std::size_t max_matrix_spectrum_bytes = 256 * 1024 * 1024) :
            bn254_small_prime_fft_backend(checked_params(max_transform_size, max_spectrum_bytes),
                                          max_matrix_spectrum_bytes) {
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

        /**
         * Multiply a row-major 2-by-2 matrix by a 2-by-Columns matrix, sharing input transforms and
         * reconstructing each sum of two products once. Columns is one for a polynomial pair.
         * Outputs must be distinct and all pointers nonnull; outputs may alias any input.
         * Return false without changing outputs when the shape, CRT bound or workspace budget does
         * not permit fusion. Small or heavily unbalanced products retain the ordinary dispatch.
         */
        template<std::size_t Columns>
            requires(Columns == 1 || Columns == 2)
        bool try_multiply_matrix(const std::array<polynomial_type *, 2 * Columns> &output,
                                 const std::array<const polynomial_type *, 4> &left,
                                 const std::array<const polynomial_type *, 2 * Columns> &right) {
            if (plans_.empty() || max_matrix_spectrum_bytes_ == 0) {
                return false;
            }
            for (std::size_t i = 0; i < output.size(); ++i) {
                if (output[i] == nullptr ||
                    std::find(output.begin(), output.begin() + i, output[i]) != output.begin() + i) {
                    return false;
                }
            }
            std::array<const polynomial_type *, 4 + 2 * Columns> operands;
            std::copy(left.begin(), left.end(), operands.begin());
            std::copy(right.begin(), right.end(), operands.begin() + 4);
            for (const auto *operand : operands) {
                if (operand == nullptr || operand->size() <= 2 || math::is_zero(*operand)) {
                    return false;
                }
            }

            std::array<std::size_t, 2 * Columns> result_sizes {};
            std::size_t max_product_size = 0;
            std::size_t min_product_size = max_transform_size_;
            const std::size_t term_limit = std::max(std::size_t(1), max_transform_size_ / 2);
            for (std::size_t row = 0; row < 2; ++row) {
                for (std::size_t column = 0; column < Columns; ++column) {
                    std::size_t terms = 0;
                    for (std::size_t inner = 0; inner < 2; ++inner) {
                        const std::size_t a = left[2 * row + inner]->size();
                        const std::size_t b = right[Columns * inner + column]->size();
                        if (a > max_transform_size_ || b > max_transform_size_ - (a - 1)) {
                            return false;
                        }
                        // Each output coefficient sums at most min(a,b) tower products from this
                        // convolution. Their combined bound must fit the existing signed CRT plan:
                        // 102 * (min(a,b) + min(c,d)) * (p-1)^2 <= its configured coefficient bound.
                        const std::size_t product_terms = std::min(a, b);
                        if (product_terms > term_limit - terms) {
                            return false;
                        }
                        terms += product_terms;
                        const std::size_t product_size = a + b - 1;
                        auto &result_size = result_sizes[row * Columns + column];
                        result_size = std::max(result_size, product_size);
                        max_product_size = std::max(max_product_size, product_size);
                        min_product_size = std::min(min_product_size, product_size);
                    }
                }
            }
            const std::size_t size = std::bit_ceil(max_product_size);
            if (std::bit_ceil(min_product_size) < size / 2) {
                return false;
            }

            // Equal input objects share storage as well as their conversion and forward transforms.
            std::array<const polynomial_type *, 4 + 2 * Columns> unique_operands {};
            std::array<std::size_t, 4 + 2 * Columns> spectrum_indices;
            std::size_t unique_count = 0;
            for (std::size_t i = 0; i < operands.size(); ++i) {
                auto end = unique_operands.begin() + unique_count;
                auto found = std::find(unique_operands.begin(), end, operands[i]);
                spectrum_indices[i] = found - unique_operands.begin();
                if (found == end) {
                    unique_operands[unique_count++] = operands[i];
                }
            }
            const std::size_t cells = field_type::arity * prime_count() * size;
            if (cells > max_matrix_spectrum_bytes_ / sizeof(double) / unique_count) {
                return false;
            }
            const std::size_t total_cells = unique_count * cells;
            if (matrix_spectra_.capacity() < total_cells) {
                // Avoid geometric vector growth beyond the configured spectrum budget.
                std::vector<double>(total_cells).swap(matrix_spectra_);
            } else {
                matrix_spectra_.resize(total_cells);
            }
            for (std::size_t i = 0; i < unique_count; ++i) {
                const auto &operand = *unique_operands[i];
                prepare_spectrum(std::span<double>(matrix_spectra_).subspan(i * cells, cells),
                                 std::span<const value_type>(operand.data(), operand.size()), size);
            }

            left_spectrum_.resize(cells);
            std::array<polynomial_type, 2 * Columns> results;
            const std::size_t plan_offset = std::countr_zero(size) * prime_count();
            for (std::size_t row = 0; row < 2; ++row) {
                for (std::size_t column = 0; column < Columns; ++column) {
                    for (std::size_t prime = 0; prime < prime_count(); ++prime) {
                        const auto &plan = plans_[plan_offset + prime];
                        const std::size_t count = field_type::arity * size;
                        const std::size_t offset = prime * count;
                        const auto spectrum = [&](std::size_t operand) {
                            return std::span<const double>(matrix_spectra_)
                                .subspan(spectrum_indices[operand] * cells + offset, count);
                        };
                        auto destination = std::span<double>(left_spectrum_).subspan(offset, count);
                        pointwise_[prime].multiply_sum(destination, spectrum(2 * row), spectrum(4 + column),
                                                       spectrum(2 * row + 1), spectrum(4 + Columns + column), size,
                                                       plan.size_inverse());
                        for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                            plan.inverse_fft(destination.subspan(coordinate * size, size), false);
                        }
                    }
                    const std::size_t index = row * Columns + column;
                    auto &result = results[index];
                    result.resize(result_sizes[index]);
                    reconstructor_.reconstruct(std::span<value_type>(result.data(), result.size()), left_spectrum_,
                                               size);
                    math::condense(result);
                }
            }
            // All inputs have been consumed and all allocations/reconstruction have succeeded.
            for (std::size_t i = 0; i < output.size(); ++i) {
                output[i]->get_storage().swap(results[i].get_storage());
            }
            return true;
        }

        /**
         * Prepare the fixed right operand for multiply_low, selecting the same prefixes and smallest
         * transform as an ordinary product. Positive precision with no variable coefficients returns
         * nullopt. Nonzero preparations also require a live backend and capacity for the complete
         * prefix product. Zero precision and zero fixed prefixes need no transform or spectrum allocation.
         * The source is a coefficient polynomial; trailing zero coefficients are permitted.
         */
        std::optional<prepared_low_product_type> prepare_low_product(const polynomial_type &fixed,
                                                                     std::size_t variable_coefficient_count,
                                                                     std::size_t coefficient_count) {
            const std::size_t prefix_size = std::min(fixed.size(), coefficient_count);
            if (coefficient_count == 0) {
                return prepared_low_product_type(coefficient_count, prefix_size);
            }
            const std::size_t variable_size = std::min(variable_coefficient_count, coefficient_count);
            if (variable_size == 0) {
                return std::nullopt;
            }
            if (math::is_zero(fixed.begin(), fixed.begin() + prefix_size)) {
                return prepared_low_product_type(coefficient_count, prefix_size);
            }
            if (plans_.empty() || variable_size > max_transform_size_ ||
                prefix_size > max_transform_size_ - (variable_size - 1)) {
                return std::nullopt;
            }
            const std::size_t size = std::bit_ceil(variable_size + prefix_size - 1);
            return prepared_low_product_type(
                coefficient_count, prefix_size, size, prime_count(),
                prepare_spectrum(std::span<const value_type>(fixed.data(), prefix_size), size));
        }

        /**
         * Compute left*fixed mod X^coefficient_count without converting or transforming fixed again.
         * Return false, leaving output unchanged, for mismatched precision, a moved-from nonzero
         * preparation, or an incompatible backend/shape. Nonzero products require the same selected
         * FFT length and catalogue prefix; zero products need neither. Successful output is canonical
         * and may alias left. The preparation is never modified or copied into working storage.
         *
         * Matching primes and a product length fitting the current backend suffice for CRT safety:
         * its existing coefficient bound covers every linear product that fits, even when the
         * preparation came from a backend with a different configured maximum.
         */
        bool try_multiply_low_prepared(polynomial_type &output, const polynomial_type &left,
                                       const prepared_low_product_type &prepared, std::size_t coefficient_count) {
            if (coefficient_count != prepared.coefficient_count_) {
                return false;
            }
            if (prepared.transform_size_ == 0) {
                output.assign(1, value_type::zero());
                return true;
            }
            if (prepared.spectrum_.empty()) {
                return false;
            }
            const std::size_t left_size = std::min(left.size(), coefficient_count);
            if (math::is_zero(left.begin(), left.begin() + left_size)) {
                output.assign(1, value_type::zero());
                return true;
            }
            if (plans_.empty() || prepared.prime_count_ != prime_count() || left_size > max_transform_size_ ||
                prepared.prefix_size_ > max_transform_size_ - (left_size - 1)) {
                return false;
            }
            const std::size_t product_size = left_size + prepared.prefix_size_ - 1;
            const std::size_t size = std::bit_ceil(product_size);
            if (size != prepared.transform_size_) {
                return false;
            }
            left_spectrum_.resize(field_type::arity * prime_count() * size);
            converter_.convert(left_spectrum_, std::span<const value_type>(left.data(), left_size), size);
            finish_product<false, true>(output, size, std::min(product_size, coefficient_count), prepared.spectrum_);
            return true;
        }

        /**
         * Prepare the complete divisor H for cyclic remainder reconstruction in the reduced-product
         * range degree(V) <= 2*n-2, n=degree(H). Cap the expected quotient length at n-1. Use the smallest
         * N >= n only when it is shorter than the supported ordinary low-product transform. Constant,
         * linear and monomial divisors, unsuitable sizes and moved-from backends return nullopt.
         *
         * Fold H in Fp12 before conversion: if N==n, add its leading coefficient into cell zero.
         * Canonical field addition keeps every converted coordinate in [0,p-1], preserving the CRT
         * bound used by cyclic_remainder_size. The preparation owns its spectrum and cached metadata.
         *
         * @pre divisor is a nonempty canonical coefficient polynomial.
         */
        std::optional<prepared_cyclic_remainder_type> prepare_cyclic_remainder(const polynomial_type &divisor,
                                                                               std::size_t quotient_coefficient_count) {
            if (divisor.size() < 3) {
                return std::nullopt;
            }
            const std::size_t n = divisor.size() - 1;
            const std::size_t size = cyclic_remainder_size(n, std::min(quotient_coefficient_count, n - 1));
            if (size == 0 || math::is_zero(divisor.begin(), divisor.begin() + n)) {
                return std::nullopt;
            }
            std::span<const value_type> coefficients(divisor.data(), divisor.size());
            if (size == n) {
                result_.get_storage().assign(divisor.begin(), divisor.begin() + n);
                result_[0] += divisor[n];
                coefficients = std::span<const value_type>(result_.data(), n);
            }
            return prepared_cyclic_remainder_type(n, size, prime_count(), prepare_spectrum(coefficients, size));
        }

        /**
         * Reconstruct R=fold_N(V)-(Q*H mod (X^N-1)). Since degree(R)<n<=N, wrapped high coefficients
         * cancel. Return false without changing output outside the reduced-product degree range,
         * without a compatible shorter transform/catalogue prefix, or for a moved-from preparation.
         * Successful output is canonical and may alias either input; the preparation is unchanged.
         *
         * @pre dividend and quotient are nonempty canonical coefficient polynomials. quotient is the
         * exact quotient of dividend by the divisor used to construct prepared.
         */
        bool try_cyclic_remainder(polynomial_type &output, const polynomial_type &dividend,
                                  const polynomial_type &quotient, const prepared_cyclic_remainder_type &prepared) {
            const std::size_t n = prepared.divisor_degree_;
            if (prepared.spectrum_.empty() || dividend.size() <= n || dividend.size() - n > n - 1 ||
                quotient.size() != dividend.size() - n || prepared.prime_count_ != prime_count()) {
                return false;
            }
            const std::size_t size = cyclic_remainder_size(n, quotient.size());
            if (size == 0 || size != prepared.transform_size_) {
                return false;
            }
            // Q has fewer than n coefficients, so it needs only zero padding, never cyclic folding.
            left_spectrum_.resize(field_type::arity * prime_count() * size);
            converter_.convert(left_spectrum_, std::span<const value_type>(quotient.data(), quotient.size()), size);
            finish_product<false, true>(output, size, size, prepared.spectrum_,
                                        std::span<const value_type>(dividend.data(), dividend.size()));
            return true;
        }

    private:
        std::size_t cyclic_remainder_size(std::size_t n, std::size_t quotient_size) const {
            if (plans_.empty() || n < 2 || quotient_size == 0 || quotient_size > n - 1 || n > max_transform_size_ ||
                quotient_size > max_transform_size_ - (n - 1)) {
                return 0;
            }
            // Require the ordinary Q*(H mod X^n) product to fit, as in the existing cyclic backend.
            // Then 2*quotient_size <= quotient_size+n-1 <= Nmax. Each cyclic coefficient contains at
            // most quotient_size products of canonical Fp12 values, so the existing Nmax/2 CRT bound
            // remains sufficient even across backend instances. Do not infer this bound from N alone.
            const std::size_t size = std::bit_ceil(n);
            return size < std::bit_ceil(quotient_size + n - 1) ? size : 0;
        }

        std::vector<double> prepare_spectrum(std::span<const value_type> coefficients, std::size_t size) const {
            // Allocate just this spectrum, so a small preparation never retains larger scratch capacity.
            std::vector<double> spectrum(field_type::arity * prime_count() * size);
            prepare_spectrum(spectrum, coefficients, size);
            return spectrum;
        }

        void prepare_spectrum(std::span<double> spectrum, std::span<const value_type> coefficients,
                              std::size_t size) const {
            converter_.convert(spectrum, coefficients, size);
            const std::size_t plan_offset = std::countr_zero(size) * prime_count();
            for (std::size_t prime = 0; prime < prime_count(); ++prime) {
                for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                    plans_[plan_offset + prime].fft(
                        spectrum.subspan((prime * field_type::arity + coordinate) * size, size));
                }
            }
        }

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

        explicit bn254_small_prime_fft_backend(const params_type &params, std::size_t max_matrix_spectrum_bytes) :
            max_transform_size_(params.transform_size), max_matrix_spectrum_bytes_(max_matrix_spectrum_bytes),
            converter_(params.crt.prime_count), reconstructor_(params.crt, params.coefficient_bound) {
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
            const std::size_t result_size = std::min(product_size, coefficient_count);
            // Euclidean quotients are commonly linear. Their matrix updates need only one or two
            // coefficient products per output, without conversion, transforms or reconstruction.
            // Keep the capacity check above: direct arithmetic has the same supported shapes.
            if (std::min(left_size, right_size) <= 2) {
                const bool left_is_shorter = left_size <= right_size;
                schoolbook_backend<value_type> {}.multiply_low(output, left_is_shorter ? left : right,
                                                               left_is_shorter ? right : left, result_size);
                return;
            }
            const std::size_t size = std::bit_ceil(product_size);
            const std::size_t channels = field_type::arity * prime_count();
            // Constructor checks guarantee that even the largest spectrum's byte count fits size_t.
            left_spectrum_.resize(channels * size);
            converter_.convert(left_spectrum_, std::span<const value_type>(left.data(), left_size), size);
            if constexpr (!Square) {
                right_spectrum_.resize(channels * size);
                converter_.convert(right_spectrum_, std::span<const value_type>(right.data(), right_size), size);
            }
            finish_product<Square>(output, size, result_size);
        }

        // Keep ordinary multiplication, squaring and prepared reuse on the same inverse/CRT path.
        // A prepared right spectrum is read directly and skips all of the right operand's forward FFTs.
        template<bool Square, bool Prepared = false>
        void finish_product(polynomial_type &output, std::size_t size, std::size_t result_size,
                            std::span<const double> prepared_spectrum = {},
                            std::span<const value_type> cyclic_dividend = {}) {
            const std::size_t plan_offset = std::countr_zero(size) * prime_count();
            for (std::size_t prime = 0; prime < prime_count(); ++prime) {
                const auto &plan = plans_[plan_offset + prime];
                const std::size_t offset = prime * field_type::arity * size;
                auto a = std::span<double>(left_spectrum_).subspan(offset, field_type::arity * size);
                std::span<const double> b;
                if constexpr (!Square) {
                    if constexpr (Prepared) {
                        b = prepared_spectrum.subspan(offset, field_type::arity * size);
                    } else {
                        b = std::span<const double>(right_spectrum_).subspan(offset, field_type::arity * size);
                    }
                }
                for (std::size_t coordinate = 0; coordinate < field_type::arity; ++coordinate) {
                    plan.fft(a.subspan(coordinate * size, size));
                    if constexpr (!Square && !Prepared) {
                        plan.fft(std::span<double>(right_spectrum_).subspan(offset + coordinate * size, size));
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
            if (!cyclic_dividend.empty()) {
                // Here result_size==size and degree(V)<=2*n-2<2*size, so at most one high coefficient
                // folds into each cell. Keep the complete difference: any uncancelled tail stays visible.
                for (std::size_t i = 0; i < size; ++i) {
                    value_type folded = i < cyclic_dividend.size() ? cyclic_dividend[i] : value_type::zero();
                    if (i + size < cyclic_dividend.size()) {
                        folded += cyclic_dividend[i + size];
                    }
                    result_[i] = folded - result_[i];
                }
            }
            math::condense(result_);
            output.get_storage().swap(result_.get_storage());
        }

        std::size_t max_transform_size_;
        std::size_t max_matrix_spectrum_bytes_;
        math::detail::bn254_small_prime_fft_converter converter_;
        math::detail::bn254_small_prime_fft_reconstructor reconstructor_;
        std::vector<math::detail::small_prime_fft_plan> plans_;
        std::vector<math::detail::bn254_small_prime_fft_pointwise> pointwise_;
        std::vector<double> left_spectrum_;
        std::vector<double> right_spectrum_;
        std::vector<double> matrix_spectra_;
        polynomial_type result_;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_BN254_SMALL_PRIME_FFT_BACKEND_HPP
