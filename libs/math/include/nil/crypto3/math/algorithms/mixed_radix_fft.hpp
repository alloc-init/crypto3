//---------------------------------------------------------------------------//
// Copyright (c) 2026 Riccardo Abbate
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

#ifndef CRYPTO3_MATH_MIXED_RADIX_FFT_HPP
#define CRYPTO3_MATH_MIXED_RADIX_FFT_HPP

#include <array>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/type_traits.hpp>

#include <nil/crypto3/math/algorithms/unity_root.hpp>

namespace nil {
    namespace crypto3 {
        namespace math {

            namespace detail {
                template<typename ValueType, typename ScalarType>
                concept MontgomeryCoordinateFieldElement =
                    algebra::FieldElementWithCoordinates<ValueType> &&
                    requires(ValueType &value, const ScalarType &scalar, std::size_t index) {
                        value.coordinate(index).data.backend().base_data();
                        scalar.data.backend().base_data();
                    };
            }    // namespace detail

            /**
             * Stores the field-only data for exact-size mixed-radix transforms.
             *
             * The plan factors the transform size into prime radices and caches
             * every forward and inverse power of a primitive root of that size.
             * It is independent of the type of values that a transform operates on.
             */
            template<typename FieldType>
            class mixed_radix_fft_plan {
                static_assert(algebra::is_field<FieldType>::value, "FieldType must be a field");

                using field_value_type = typename FieldType::value_type;

            public:
                explicit mixed_radix_fft_plan(std::size_t size) :
                    size_(validate_size(size)), radices_(detail::prime_factors(size_)),
                    omega_(unity_root<FieldType>(size_)), size_inverse_(field_value_type(size_).inversed()),
                    forward_powers_(root_powers(omega_, size_)),
                    inverse_powers_(root_powers(forward_powers_[size_ - 1], size_)) {
                }

                std::size_t size() const {
                    return size_;
                }

                const std::vector<std::size_t> &radices() const {
                    return radices_;
                }

                const field_value_type &omega() const {
                    return omega_;
                }

                /**
                 * Replace coefficients with their evaluations at successive powers of omega.
                 * Short inputs are zero-padded to the plan size; oversized inputs are rejected.
                 */
                template<typename ValueType>
                void fft(std::vector<ValueType> &values) const {
                    if (values.size() > size_) {
                        throw std::invalid_argument("mixed_radix_fft_plan::fft: input exceeds plan size");
                    }

                    values.resize(size_, ValueType::zero());
                    std::vector<ValueType> scratch(size_, ValueType::zero());
                    transform_recursive(values, 0, 1, scratch, 0, size_, 0, 1, forward_powers_);
                    values.swap(scratch);
                }

                /**
                 * Recover coefficients from evaluations at successive powers of omega.
                 * Short inputs are zero-padded to the plan size; oversized inputs are rejected.
                 */
                template<typename ValueType>
                void inverse_fft(std::vector<ValueType> &values) const {
                    if (values.size() > size_) {
                        throw std::invalid_argument("mixed_radix_fft_plan::inverse_fft: input exceeds plan size");
                    }

                    values.resize(size_, ValueType::zero());
                    std::vector<ValueType> scratch(size_, ValueType::zero());
                    transform_recursive(values, 0, 1, scratch, 0, size_, 0, 1, inverse_powers_);
                    for (ValueType &value : scratch) {
                        value = value * size_inverse_;
                    }
                    values.swap(scratch);
                }

            private:
                // Reject zero before the size is used to construct the rest of the plan.
                static std::size_t validate_size(std::size_t size) {
                    if (size == 0) {
                        throw std::invalid_argument("mixed_radix_fft_plan: expected size > 0");
                    }
                    return size;
                }

                // Build [1, root, root^2, ..., root^(size - 1)] using one multiplication per new power.
                static std::vector<field_value_type> root_powers(const field_value_type &root, std::size_t size) {
                    std::vector<field_value_type> powers(size, field_value_type::one());
                    for (std::size_t i = 1; i < size; ++i) {
                        powers[i] = powers[i - 1] * root;
                    }
                    return powers;
                }

                // Apply the backend inner product independently to each base-field coordinate. The products remain
                // unreduced until the complete radix sum has been accumulated, replacing one Montgomery reduction per
                // term with one reduction per output coordinate.
                template<typename ValueType>
                    requires detail::MontgomeryCoordinateFieldElement<ValueType, field_value_type>
                void lazy_montgomery_radix_dft(std::vector<ValueType> &workspace, std::size_t workspace_offset,
                                               std::size_t workspace_stride, std::vector<ValueType> &output,
                                               std::size_t output_offset, std::size_t output_stride, std::size_t radix,
                                               std::size_t root_stride,
                                               const std::vector<field_value_type> &powers) const {
                    using base_value_type =
                        std::remove_cvref_t<decltype(std::declval<ValueType &>().coordinate(std::size_t(0)))>;
                    using value_field_type = typename ValueType::field_type;
                    constexpr std::size_t coordinate_count = value_field_type::arity;
                    using backend_type =
                        std::remove_reference_t<decltype(std::declval<base_value_type &>().data.backend().base_data())>;
                    using scalar_backend_type = std::remove_reference_t<
                        decltype(std::declval<field_value_type &>().data.backend().base_data())>;
                    static_assert(std::is_same_v<backend_type, scalar_backend_type>,
                                  "FFT values and roots must use the same base-field backend");

                    std::array<std::vector<backend_type>, coordinate_count> coordinate_inputs;
                    for (std::vector<backend_type> &inputs : coordinate_inputs) {
                        inputs.resize(radix);
                    }
                    for (std::size_t input_index = 0; input_index < radix; ++input_index) {
                        ValueType &value = workspace[workspace_offset + input_index * workspace_stride];
                        for (std::size_t coordinate_index = 0; coordinate_index < coordinate_count;
                             ++coordinate_index) {
                            coordinate_inputs[coordinate_index][input_index] =
                                value.coordinate(coordinate_index).data.backend().base_data();
                        }
                    }

                    std::vector<backend_type> scalar_inputs(radix);
                    for (std::size_t output_index = 0; output_index < radix; ++output_index) {
                        std::size_t power_index = 0;
                        const std::size_t power_step = root_stride * output_index;
                        for (std::size_t input_index = 0; input_index < radix; ++input_index) {
                            scalar_inputs[input_index] = powers[power_index].data.backend().base_data();
                            power_index += power_step;
                            if (power_index >= size_) {
                                power_index -= size_;
                            }
                        }

                        ValueType sum = ValueType::zero();
                        for (std::size_t coordinate_index = 0; coordinate_index < coordinate_count;
                             ++coordinate_index) {
                            backend_type coordinate_sum;
                            base_value_type::modulus_params.get_mod_obj().montgomery_inner_product(
                                coordinate_sum, coordinate_inputs[coordinate_index].begin(),
                                coordinate_inputs[coordinate_index].end(), scalar_inputs.begin());
                            sum.coordinate(coordinate_index).data.backend().base_data() = coordinate_sum;
                        }
                        output[output_offset + output_index * output_stride] = sum;
                    }
                }

                // Evaluate one prime-radix DFT from temporary input slots into strided output slots.
                // Direct radix kernels make the full transform O(size * sum(radices)); this is efficient for
                // small and moderate prime factors, but sizes with very large prime factors need another algorithm.
                template<typename ValueType>
                void radix_dft(std::vector<ValueType> &workspace, std::size_t workspace_offset,
                               std::size_t workspace_stride, std::vector<ValueType> &output, std::size_t output_offset,
                               std::size_t output_stride, std::size_t radix, std::size_t root_stride,
                               const std::vector<field_value_type> &powers) const {
                    // Constructing backend input arrays is more expensive than eager reduction for the small radix-2
                    // and radix-3 kernels. Larger kernels amortize that setup over enough products to benefit from
                    // accumulating each base-field coordinate before a single Montgomery reduction.
                    constexpr std::size_t lazy_montgomery_radix_threshold = 8;
                    if (radix >= lazy_montgomery_radix_threshold) {
                        if constexpr (detail::MontgomeryCoordinateFieldElement<ValueType, field_value_type>) {
                            lazy_montgomery_radix_dft(workspace, workspace_offset, workspace_stride, output,
                                                      output_offset, output_stride, radix, root_stride, powers);
                            return;
                        }
                    }

                    // Generic values and small radices use the ordinary field operations.
                    for (std::size_t output_index = 0; output_index < radix; ++output_index) {
                        ValueType sum = ValueType::zero();
                        std::size_t power_index = 0;
                        const std::size_t power_step = root_stride * output_index;

                        for (std::size_t input_index = 0; input_index < radix; ++input_index) {
                            sum += workspace[workspace_offset + input_index * workspace_stride] * powers[power_index];
                            power_index += power_step;
                            if (power_index >= size_) {
                                power_index -= size_;
                            }
                        }

                        output[output_offset + output_index * output_stride] = sum;
                    }
                }

                // Recursively transform residue classes modulo the next radix, then combine them in standard order.
                template<typename ValueType>
                void transform_recursive(std::vector<ValueType> &input, std::size_t input_offset,
                                         std::size_t input_stride, std::vector<ValueType> &output,
                                         std::size_t output_offset, std::size_t transform_size, std::size_t radix_index,
                                         std::size_t root_stride, const std::vector<field_value_type> &powers) const {
                    if (transform_size == 1) {
                        output[output_offset] = input[input_offset];
                        return;
                    }

                    const std::size_t radix = radices_[radix_index];
                    const std::size_t subtransform_size = transform_size / radix;

                    for (std::size_t residue = 0; residue < radix; ++residue) {
                        transform_recursive(input, input_offset + residue * input_stride, input_stride * radix, output,
                                            output_offset + residue * subtransform_size, subtransform_size,
                                            radix_index + 1, root_stride * radix, powers);
                    }

                    for (std::size_t subtransform_index = 0; subtransform_index < subtransform_size;
                         ++subtransform_index) {
                        std::size_t power_index = 0;
                        const std::size_t power_step = root_stride * subtransform_index;

                        for (std::size_t residue = 0; residue < radix; ++residue) {
                            input[input_offset + residue * input_stride] =
                                output[output_offset + residue * subtransform_size + subtransform_index] *
                                powers[power_index];
                            power_index += power_step;
                            if (power_index >= size_) {
                                power_index -= size_;
                            }
                        }

                        radix_dft(input, input_offset, input_stride, output, output_offset + subtransform_index,
                                  subtransform_size, radix, root_stride * subtransform_size, powers);
                    }
                }

                std::size_t size_;
                std::vector<std::size_t> radices_;
                field_value_type omega_;
                field_value_type size_inverse_;
                std::vector<field_value_type> forward_powers_;
                std::vector<field_value_type> inverse_powers_;
            };

        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MATH_MIXED_RADIX_FFT_HPP
