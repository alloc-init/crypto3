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

#ifndef CRYPTO3_MATH_MIXED_RADIX_BACKEND_HPP
#define CRYPTO3_MATH_MIXED_RADIX_BACKEND_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/math/algorithms/mixed_radix_fft.hpp>
#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>
#include <nil/crypto3/math/polynomial/backends/polynomial_backend.hpp>

namespace nil::crypto3::math::polynomial_arithmetic {

    /**
     * Mixed-radix multiplication backend with an explicitly configured maximum transform order.
     * The backend caches plans for the divisors of that order and uses the smallest plan that
     * contains each product. It reuses its workspace across operations and is therefore intended
     * for sequential use. Concurrent operations require separate backend instances.
     *
     * The configured transform must contain the complete product of the coefficients
     * used by an operation. multiply_low discards input coefficients that cannot
     * affect the requested low coefficients before computing that product.
     */
    template<typename RootFieldType, typename ValueType = typename RootFieldType::value_type>
    class mixed_radix_backend {
    public:
        using value_type = ValueType;
        using polynomial_type = math::polynomial<value_type>;

        /**
         * Owned, immutable preparation of a fixed operand for multiply_low. The spectrum contains the FFT of
         * the first min(fixed.size(), coefficient_count) coefficients, with enough zero padding for a linear
         * product. It does not refer to the source polynomial, backend, FFT plan, or scratch storage.
         */
        class prepared_low_product_type {
        public:
            std::size_t coefficient_count() const {
                return coefficient_count_;
            }

            // Zero means the requested fixed prefix is zero, so no transform was needed.
            std::size_t transform_size() const {
                return spectrum_.size();
            }

            std::size_t storage_bytes() const {
                return spectrum_.capacity() * sizeof(value_type);
            }

        private:
            friend class mixed_radix_backend;

            prepared_low_product_type(std::size_t coefficient_count, std::size_t prefix_size, bool zero,
                                      std::vector<value_type> spectrum = {}) :
                coefficient_count_(coefficient_count), prefix_size_(prefix_size), zero_(zero),
                spectrum_(std::move(spectrum)) {
            }

            std::size_t coefficient_count_;
            std::size_t prefix_size_;
            bool zero_;
            std::vector<value_type> spectrum_;
        };

        explicit mixed_radix_backend(std::size_t transform_size) {
            // Every divisor of a valid transform order also has a root of unity. Caching those plans lets algorithms
            // whose operands grow in stages use a proportionally sized transform at each stage instead of repeatedly
            // paying for the configured maximum transform.
            const std::vector<std::size_t> sizes = divisor_transform_sizes(transform_size);
            plans_.reserve(sizes.size());
            for (const std::size_t size : sizes) {
                plans_.emplace_back(size);
            }
        }

        void multiply(polynomial_type &output, const polynomial_type &left, const polynomial_type &right) {
            if (math::is_zero(left) || math::is_zero(right)) {
                set_zero(output);
                return;
            }

            const std::size_t result_size = left.size() + right.size() - 1;
            multiply_prefixes(output, left, left.size(), right, right.size(), result_size);
        }

        void square(polynomial_type &output, const polynomial_type &input) {
            if (math::is_zero(input)) {
                set_zero(output);
                return;
            }

            const std::size_t result_size = 2 * input.size() - 1;
            const mixed_radix_fft_plan<RootFieldType> &plan = plan_for(result_size);

            polynomial_type transformed = input;
            plan.fft(transformed.get_storage(), workspace_);
            for (value_type &value : transformed) {
                value = value * value;
            }

            finish_transform(output, transformed, result_size, plan);
        }

        void multiply_low(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                          std::size_t coefficient_count) {
            if (coefficient_count == 0) {
                set_zero(output);
                return;
            }

            const std::size_t left_size = std::min(left.size(), coefficient_count);
            const std::size_t right_size = std::min(right.size(), coefficient_count);
            if (math::is_zero(left.begin(), left.begin() + left_size) ||
                math::is_zero(right.begin(), right.begin() + right_size)) {
                set_zero(output);
                return;
            }

            const std::size_t prefix_product_size = left_size + right_size - 1;
            const std::size_t result_size = std::min(coefficient_count, prefix_product_size);
            multiply_prefixes(output, left, left_size, right, right_size, result_size);
        }

        /**
         * Prepare a fixed right operand for a low product with a variable operand of variable_coefficient_count
         * coefficients. This uses the same prefixes and smallest fitting transform as multiply_low.
         * Returns nullopt if a positive precision has no variable coefficients or no plan can hold the complete
         * prefix product. Zero precision and a zero fixed prefix need no FFT and produce a zero preparation.
         *
         * @pre fixed is a nonempty canonical coefficient polynomial.
         */
        std::optional<prepared_low_product_type> prepare_low_product(const polynomial_type &fixed,
                                                                     std::size_t variable_coefficient_count,
                                                                     std::size_t coefficient_count) {
            const std::size_t prefix_size = std::min(fixed.size(), coefficient_count);
            if (coefficient_count == 0) {
                return prepared_low_product_type(coefficient_count, prefix_size, true);
            }
            const std::size_t variable_size = std::min(variable_coefficient_count, coefficient_count);
            if (variable_size == 0) {
                return std::nullopt;
            }
            if (math::is_zero(fixed.begin(), fixed.begin() + prefix_size)) {
                return prepared_low_product_type(coefficient_count, prefix_size, true);
            }
            if (variable_size > std::numeric_limits<std::size_t>::max() - (prefix_size - 1)) {
                return std::nullopt;
            }
            const auto *plan = find_plan(variable_size + prefix_size - 1);
            if (!plan) {
                return std::nullopt;
            }

            std::vector<value_type> spectrum(fixed.begin(), fixed.begin() + prefix_size);
            plan->fft(spectrum, workspace_);
            // FFT exchanges buffers with workspace_. A preceding larger transform can therefore leave this
            // smaller spectrum with a large allocation. Detach that excess capacity once during preparation.
            if (spectrum.capacity() > spectrum.size()) {
                std::vector<value_type>(spectrum.begin(), spectrum.end()).swap(spectrum);
            }
            return prepared_low_product_type(coefficient_count, prefix_size, false, std::move(spectrum));
        }

        /**
         * Compute left * fixed mod X^coefficient_count using an existing preparation. Returns false, leaving
         * output unchanged, if the precision differs or the current backend would select a different FFT length.
         * Different variable-operand lengths and backend instances are compatible when that selected length
         * agrees. Zero products need no matching FFT plan. Successful output is canonical and may alias left;
         * the preparation is never modified.
         *
         * @pre left is a nonempty coefficient polynomial. Trailing zero coefficients are permitted, as in the
         *      fixed-length reversal of a dividend prefix during division.
         */
        bool try_multiply_low_prepared(polynomial_type &output, const polynomial_type &left,
                                       const prepared_low_product_type &prepared, std::size_t coefficient_count) {
            if (coefficient_count != prepared.coefficient_count_) {
                return false;
            }
            if (prepared.zero_) {
                set_zero(output);
                return true;
            }
            // A moved-from nonzero preparation no longer owns a spectrum.
            if (prepared.spectrum_.empty()) {
                return false;
            }
            const std::size_t left_size = std::min(left.size(), coefficient_count);
            if (math::is_zero(left.begin(), left.begin() + left_size)) {
                set_zero(output);
                return true;
            }
            if (left_size > std::numeric_limits<std::size_t>::max() - (prepared.prefix_size_ - 1)) {
                return false;
            }
            const std::size_t product_size = left_size + prepared.prefix_size_ - 1;
            const auto *plan = find_plan(product_size);
            if (!plan || plan->size() != prepared.transform_size()) {
                return false;
            }

            polynomial_type transformed(left.begin(), left.begin() + left_size);
            plan->fft(transformed.get_storage(), workspace_);
            finish_product(output, transformed, prepared.spectrum_, std::min(coefficient_count, product_size), *plan);
            return true;
        }

    private:
        static std::vector<std::size_t> divisor_transform_sizes(std::size_t transform_size) {
            if (transform_size == 0) {
                throw std::invalid_argument("mixed_radix_backend: expected transform size > 0");
            }

            const std::vector<std::size_t> factors = math::detail::prime_factors(transform_size);
            std::vector<std::size_t> sizes = {1};
            for (std::size_t factor_index = 0; factor_index < factors.size();) {
                // If the next prime occurs e times, combine each divisor built so far with p, ..., p^e.
                const std::size_t factor = factors[factor_index];
                std::size_t next_factor_index = factor_index;
                while (next_factor_index < factors.size() && factors[next_factor_index] == factor) {
                    ++next_factor_index;
                }

                const std::size_t existing_size_count = sizes.size();
                std::size_t factor_power = 1;
                for (std::size_t exponent = factor_index; exponent < next_factor_index; ++exponent) {
                    factor_power *= factor;
                    for (std::size_t i = 0; i < existing_size_count; ++i) {
                        sizes.push_back(sizes[i] * factor_power);
                    }
                }
                factor_index = next_factor_index;
            }

            std::sort(sizes.begin(), sizes.end());
            return sizes;
        }

        static void set_zero(polynomial_type &output) {
            output.assign(1, value_type::zero());
        }

        const mixed_radix_fft_plan<RootFieldType> *find_plan(std::size_t result_size) const {
            // Plans are sorted by order, so the first plan large enough for the product introduces the least padding.
            const auto plan = std::lower_bound(plans_.begin(), plans_.end(), result_size,
                                               [](const mixed_radix_fft_plan<RootFieldType> &candidate,
                                                  std::size_t size) { return candidate.size() < size; });
            return plan == plans_.end() ? nullptr : &*plan;
        }

        const mixed_radix_fft_plan<RootFieldType> &plan_for(std::size_t result_size) const {
            const auto *plan = find_plan(result_size);
            if (!plan) {
                throw std::invalid_argument("mixed_radix_backend: transform is too small for the product");
            }
            return *plan;
        }

        void multiply_prefixes(polynomial_type &output, const polynomial_type &left, std::size_t left_size,
                               const polynomial_type &right, std::size_t right_size, std::size_t result_size) {
            const mixed_radix_fft_plan<RootFieldType> &plan = plan_for(left_size + right_size - 1);

            polynomial_type transformed_left(left.begin(), left.begin() + left_size);
            polynomial_type transformed_right(right.begin(), right.begin() + right_size);
            plan.fft(transformed_left.get_storage(), workspace_);
            plan.fft(transformed_right.get_storage(), workspace_);

            finish_product(output, transformed_left, transformed_right.get_storage(), result_size, plan);
        }

        void finish_product(polynomial_type &output, polynomial_type &transformed_left,
                            const std::vector<value_type> &transformed_right, std::size_t result_size,
                            const mixed_radix_fft_plan<RootFieldType> &plan) {
            for (std::size_t i = 0; i < plan.size(); ++i) {
                transformed_left[i] = transformed_left[i] * transformed_right[i];
            }

            finish_transform(output, transformed_left, result_size, plan);
        }

        void finish_transform(polynomial_type &output, polynomial_type &transformed, std::size_t result_size,
                              const mixed_radix_fft_plan<RootFieldType> &plan) {
            plan.inverse_fft(transformed.get_storage(), workspace_);
            math::truncate(transformed, result_size);
            output = std::move(transformed);
        }

        std::vector<mixed_radix_fft_plan<RootFieldType>> plans_;
        std::vector<value_type> workspace_;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_MIXED_RADIX_BACKEND_HPP
