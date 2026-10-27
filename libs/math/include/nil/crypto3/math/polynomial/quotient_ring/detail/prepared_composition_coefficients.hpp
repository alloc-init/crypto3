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

#ifndef CRYPTO3_MATH_DETAIL_PREPARED_COMPOSITION_COEFFICIENTS_HPP
#define CRYPTO3_MATH_DETAIL_PREPARED_COMPOSITION_COEFFICIENTS_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/fp12_fast.hpp>
#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>

namespace nil::crypto3::math::detail {

    template<typename Value, typename = void>
    class prepared_composition_coefficients {
    public:
        static constexpr bool supported = false;
    };

    template<typename Value>
    class prepared_composition_coefficients<Value,
                                            std::void_t<typename Value::field_type::extension_policy::fast_params>> {
        using params = typename Value::field_type::extension_policy::fast_params;
        using kernel = algebra::fields::detail::fp12_fast::fp12_fast<params>;
        using prepared_type = typename kernel::prepared_product_operand;
        using simd_kernel = typename kernel::rank23;
        using limb = algebra::fields::detail::fp12_fast::limb;

    public:
        static constexpr bool supported = true;
        static constexpr std::size_t block_tile_size = 4;
        static constexpr std::size_t coefficient_tile_size = 4;

        // Scratch belongs to a composition call, so concurrent readers never mutate the shared cache.
        struct workspace {
            std::array<std::vector<prepared_type>, block_tile_size> scalars;
            std::array<std::vector<std::size_t>, block_tile_size> exponents;
            std::array<std::vector<limb>, block_tile_size> simd_scalars;
        };

        prepared_composition_coefficients() = default;

        template<typename Polynomial>
        prepared_composition_coefficients(const std::vector<Polynomial> &powers, std::size_t byte_limit,
                                          std::size_t simd_lanes = 8) {
            if (simd_lanes != 0 && simd_lanes != 4 && simd_lanes != 8) {
                throw std::invalid_argument("composition SIMD lanes must be 0, 4 or 8");
            }
            if (byte_limit == 0 || powers.size() < 2) {
                return;
            }
            const std::size_t terms = powers.size() - 1;    // The constant power is handled by direct addition.
            std::size_t width = 0;
            for (std::size_t exponent = 1; exponent < powers.size(); ++exponent) {
                width = std::max(width, powers[exponent].size());
            }
#if defined(CRYPTO3_HAS_AVX512IFMA)
            if constexpr (simd_kernel::supported) {
                // Cache dimensions rule out short dots and narrow polynomials. use_for_block also checks
                // the actual nonzero outer terms before selecting either prepared kernel.
                if (simd_lanes != 0 && width >= 4 && terms >= minimum_prepared_terms &&
                    prepare_simd(powers, width, terms, byte_limit, simd_lanes)) {
                    return;
                }
            }
#endif
            const std::size_t capacity = std::min(byte_limit / sizeof(prepared_type), coefficients_.max_size());
            // Division before multiplication checks both the budget and size_t overflow without allocating.
            if (width == 0 || width > capacity / terms) {
                return;
            }
            coefficients_.resize(width * terms);
            for (std::size_t coefficient = 0; coefficient < width; ++coefficient) {
                for (std::size_t exponent = 1; exponent < powers.size(); ++exponent) {
                    if (coefficient < powers[exponent].size()) {
                        coefficients_[coefficient * terms + exponent - 1] =
                            kernel::prepare_product_operand(powers[exponent][coefficient]);
                    }
                    // Default-constructed backends are zero, padding shorter powers without borrowed views.
                }
            }
            terms_ = terms;
            width_ = width;
        }

        bool available() const {
            return !coefficients_.empty() || !simd_coefficients_.empty();
        }

        std::size_t bytes() const {
            return coefficients_.size() * sizeof(prepared_type) + simd_coefficients_.size() * sizeof(limb);
        }

        std::size_t simd_lanes() const {
            return simd_lanes_;
        }

        // Short or sparse sums reuse ordinary Fp12 multiplication instead of paying for operand preparation
        // and recombination. The constant coefficient needs only an addition and does not count as a term.
        template<typename Polynomial>
        bool use_for_block(const Polynomial &outer, std::size_t block_index) const {
            if (!available()) {
                return false;
            }
            const std::size_t block_size = terms_ + 1;
            const std::size_t first = block_index * block_size;
            assert(first < outer.size());
            const std::size_t count = std::min(block_size, outer.size() - first);
            if (count <= minimum_prepared_terms) {
                return false;
            }
            std::size_t active_terms = 0;
            for (std::size_t exponent = 1; exponent < count; ++exponent) {
                if (!outer[first + exponent].is_zero() && ++active_terms == minimum_prepared_terms) {
                    return true;
                }
            }
            return false;
        }

        // Evaluate up to four adjacent outer blocks. Each scalar is prepared once and reused across all
        // output coefficients. The coefficient-major cache keeps each dot contiguous; tiles reuse the same
        // cached columns across adjacent blocks without storing a matrix of every block's output.
        template<typename Polynomial>
        void evaluate_blocks(std::array<Polynomial, block_tile_size> &output, const Polynomial &outer,
                             std::size_t first_block, std::size_t block_count, workspace &scratch) const {
            assert(available() && block_count > 0 && block_count <= block_tile_size);
#if defined(CRYPTO3_HAS_AVX512IFMA)
            if constexpr (simd_kernel::supported) {
                if (simd_lanes_ == 4) {
                    evaluate_simd_blocks<4>(output, outer, first_block, block_count, scratch);
                    return;
                }
                if (simd_lanes_ == 8) {
                    evaluate_simd_blocks<8>(output, outer, first_block, block_count, scratch);
                    return;
                }
            }
#endif
            const std::size_t block_size = terms_ + 1;
            for (std::size_t block = 0; block < block_count; ++block) {
                const std::size_t first = (first_block + block) * block_size;
                assert(first < outer.size());
                const std::size_t count = std::min(block_size, outer.size() - first);
                auto &scalars = scratch.scalars[block];
                auto &exponents = scratch.exponents[block];
                scalars.clear();
                exponents.clear();
                scalars.reserve(count - 1);
                exponents.reserve(count - 1);
                for (std::size_t exponent = 1; exponent < count; ++exponent) {
                    if (!outer[first + exponent].is_zero()) {
                        scalars.push_back(kernel::prepare_product_operand(outer[first + exponent]));
                        exponents.push_back(exponent - 1);
                    }
                }
                output[block].assign(scalars.empty() ? 1 : width_, Value::zero());
                output[block][0] = outer[first];
            }
            for (std::size_t first_coefficient = 0; first_coefficient < width_;) {
                const std::size_t end = first_coefficient + std::min(coefficient_tile_size, width_ - first_coefficient);
                for (std::size_t coefficient = first_coefficient; coefficient < end; ++coefficient) {
                    const auto select_power = [this, coefficient](std::size_t exponent) -> const prepared_type & {
                        return coefficients_[coefficient * terms_ + exponent];
                    };
                    for (std::size_t block = 0; block < block_count; ++block) {
                        const auto &scalars = scratch.scalars[block];
                        if (!scalars.empty()) {
                            const auto powers = scratch.exponents[block] | std::views::transform(select_power);
                            output[block][coefficient] +=
                                kernel::template inner_product<Value>(scalars.begin(), scalars.end(), powers.begin());
                        }
                    }
                }
                first_coefficient = end;
            }
            for (std::size_t block = 0; block < block_count; ++block) {
                condense(output[block]);
            }
        }

    private:
        // Conservative dispatch heuristic shared with the SIMD cache cutoff, not an arithmetic bound.
        static constexpr std::size_t minimum_prepared_terms = 8;
#if defined(CRYPTO3_HAS_AVX512IFMA)
        template<typename Polynomial>
        bool prepare_simd(const std::vector<Polynomial> &powers, std::size_t width, std::size_t terms,
                          std::size_t byte_limit, std::size_t lanes) {
            using radix = algebra::fields::detail::fp12_fast::montgomery_radix52<params>;
            constexpr std::size_t coordinates = simd_kernel::product_coordinate_count;
            const std::size_t panel_count = 1 + (width - 1) / lanes;
            const std::size_t stride = coordinates * 5 * lanes;
            const std::size_t capacity = std::min(byte_limit / sizeof(limb), simd_coefficients_.max_size());
            if (terms > capacity / stride || panel_count > capacity / stride / terms) {
                return false;
            }
            simd_coefficients_.resize(panel_count * stride * terms);
            // Owned [panel][coordinate][term][digit][lane] storage. Only this representation is retained;
            // 40 bytes per form replace the scalar backend object, and padding is included in the budget.
            for (std::size_t coefficient = 0; coefficient < width; ++coefficient) {
                const std::size_t panel = coefficient / lanes, lane = coefficient % lanes;
                for (std::size_t exponent = 1; exponent < powers.size(); ++exponent) {
                    if (coefficient >= powers[exponent].size()) {
                        continue;
                    }
                    const auto prepared = simd_kernel::prepare_product_operand(powers[exponent][coefficient]);
                    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
                        const auto packed = radix::from_montgomery64(prepared[coordinate]);
                        const std::size_t offset = ((panel * coordinates + coordinate) * terms + exponent - 1) * 5;
                        for (std::size_t digit = 0; digit < 5; ++digit) {
                            simd_coefficients_[(offset + digit) * lanes + lane] = packed[digit];
                        }
                    }
                }
            }
            terms_ = terms;
            width_ = width;
            simd_lanes_ = lanes;
            return true;
        }

        template<std::size_t Lanes, typename Polynomial>
        void evaluate_simd_blocks(std::array<Polynomial, block_tile_size> &output, const Polynomial &outer,
                                  std::size_t first_block, std::size_t block_count, workspace &scratch) const {
            namespace fast = algebra::fields::detail::fp12_fast;
            using radix = fast::montgomery_radix52<params>;
            constexpr std::size_t coordinates = simd_kernel::product_coordinate_count;
            bool any_active = false;
            for (std::size_t block = 0; block < block_count; ++block) {
                const std::size_t first = (first_block + block) * (terms_ + 1);
                assert(first < outer.size());
                const std::size_t count = std::min(terms_ + 1, outer.size() - first);
                auto &indices = scratch.exponents[block];
                indices.clear();
                for (std::size_t exponent = 1; exponent < count; ++exponent) {
                    if (!outer[first + exponent].is_zero()) {
                        indices.push_back(exponent - 1);
                    }
                }
                auto &scalars = scratch.simd_scalars[block];
                if (indices.size() > scalars.max_size() / (coordinates * 5)) {
                    throw std::length_error("composition SIMD scratch exceeds maximum size");
                }
                scalars.resize(indices.size() * coordinates * 5);
                for (std::size_t term = 0; term < indices.size(); ++term) {
                    const auto prepared = simd_kernel::prepare_product_operand(outer[first + indices[term] + 1]);
                    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
                        const auto packed = radix::from_montgomery64(prepared[coordinate]);
                        std::copy(packed.begin(), packed.end(),
                                  scalars.begin() + (coordinate * indices.size() + term) * 5);
                    }
                }
                any_active = any_active || !indices.empty();
                output[block].assign(indices.empty() ? 1 : width_, Value::zero());
                output[block][0] = outer[first];
            }
            for (std::size_t first = 0; any_active && first < width_;) {
                const std::size_t active_lanes = std::min(Lanes, width_ - first);
                // Keep the dot products in canonical radix-52 form until tower recombination.
                std::array<std::array<limb, coordinates * 5 * Lanes>, block_tile_size> products;
                for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
                    const limb *right =
                        simd_coefficients_.data() + ((first / Lanes * coordinates + coordinate) * terms_) * 5 * Lanes;
                    // Reuse this coordinate panel across neighboring blocks before moving to the next one.
                    for (std::size_t block = 0; block < block_count; ++block) {
                        const auto &indices = scratch.exponents[block];
                        if (indices.empty()) {
                            continue;
                        }
                        fast::montgomery_inner_product52<params, fast::ifma52_operations<Lanes>>(
                            products[block].data() + coordinate * 5 * Lanes,
                            scratch.simd_scalars[block].data() + coordinate * indices.size() * 5, right, indices.data(),
                            indices.size());
                    }
                }
                for (std::size_t block = 0; block < block_count; ++block) {
                    if (!scratch.exponents[block].empty()) {
                        std::array<limb, 12 * 5 * Lanes> combined;
                        simd_kernel::template recombine_products52<fast::ifma52_operations<Lanes>>(
                            combined.data(), products[block].data());
                        for (std::size_t lane = 0; lane < active_lanes; ++lane) {
                            Value value;
                            for (std::size_t coordinate = 0; coordinate < 12; ++coordinate) {
                                typename radix::value_type packed;
                                for (std::size_t digit = 0; digit < 5; ++digit) {
                                    packed[digit] = combined[(coordinate * 5 + digit) * Lanes + lane];
                                }
                                radix::to_montgomery64(value.coordinate(coordinate).data.backend().base_data(), packed);
                            }
                            output[block][first + lane] += value;
                        }
                    }
                }
                first += active_lanes;
            }
            for (std::size_t block = 0; block < block_count; ++block) {
                condense(output[block]);
            }
        }
#endif
        std::size_t terms_ = 0, width_ = 0, simd_lanes_ = 0;
        std::vector<prepared_type> coefficients_;
        std::vector<limb> simd_coefficients_;
    };

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_DETAIL_PREPARED_COMPOSITION_COEFFICIENTS_HPP
