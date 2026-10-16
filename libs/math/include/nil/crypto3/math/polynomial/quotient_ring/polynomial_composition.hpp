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

#ifndef CRYPTO3_MATH_POLYNOMIAL_COMPOSITION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_COMPOSITION_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <nil/crypto3/algebra/fields/detail/element/fp12_fast/fp12_fast.hpp>
#include <nil/crypto3/math/detail/integer_sqrt.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_modular_arithmetic.hpp>

namespace nil::crypto3::math {

    namespace detail {
        template<typename Value, typename = void>
        class prepared_composition_coefficients {
        public:
            static constexpr bool supported = false;
        };

        template<typename Value>
        class prepared_composition_coefficients<
            Value, std::void_t<typename Value::field_type::extension_policy::fast_params>> {
            using params = typename Value::field_type::extension_policy::fast_params;
            using kernel = algebra::fields::detail::fp12_fast::fp12_fast<params>;
            using prepared_type = typename kernel::prepared_product_operand;
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
                if constexpr (params::base_value_limb_count == 4 && params::base_field_type::modulus_bits <= 254) {
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
                if constexpr (params::base_value_limb_count == 4 && params::base_field_type::modulus_bits <= 254) {
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
                    const std::size_t end =
                        first_coefficient + std::min(coefficient_tile_size, width_ - first_coefficient);
                    for (std::size_t coefficient = first_coefficient; coefficient < end; ++coefficient) {
                        const auto select_power = [this, coefficient](std::size_t exponent) -> const prepared_type & {
                            return coefficients_[coefficient * terms_ + exponent];
                        };
                        for (std::size_t block = 0; block < block_count; ++block) {
                            const auto &scalars = scratch.scalars[block];
                            if (!scalars.empty()) {
                                const auto powers = scratch.exponents[block] | std::views::transform(select_power);
                                output[block][coefficient] += kernel::template inner_product<Value>(
                                    scalars.begin(), scalars.end(), powers.begin());
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
                constexpr std::size_t coordinates = kernel::product_coordinate_count;
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
                        const auto prepared = kernel::prepare_product_operand(powers[exponent][coefficient]);
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
                constexpr std::size_t coordinates = kernel::product_coordinate_count;
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
                        const auto prepared = kernel::prepare_product_operand(outer[first + indices[term] + 1]);
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
                    std::array<std::array<std::array<typename kernel::base_value_type, coordinates>, Lanes>,
                               block_tile_size>
                        products;
                    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
                        const limb *right = simd_coefficients_.data() +
                                            ((first / Lanes * coordinates + coordinate) * terms_) * 5 * Lanes;
                        // Reuse this coordinate panel across neighboring blocks before moving to the next one.
                        for (std::size_t block = 0; block < block_count; ++block) {
                            const auto &indices = scratch.exponents[block];
                            if (indices.empty()) {
                                continue;
                            }
                            std::array<limb, 5 * Lanes> reduced;
                            fast::montgomery_inner_product52<params, fast::ifma52_operations<Lanes>>(
                                reduced.data(), scratch.simd_scalars[block].data() + coordinate * indices.size() * 5,
                                right, indices.data(), indices.size());
                            for (std::size_t lane = 0; lane < active_lanes; ++lane) {
                                typename radix::value_type packed;
                                for (std::size_t digit = 0; digit < 5; ++digit) {
                                    packed[digit] = reduced[digit * Lanes + lane];
                                }
                                radix::to_montgomery64(products[block][lane][coordinate].data.backend().base_data(),
                                                       packed);
                            }
                        }
                    }
                    for (std::size_t block = 0; block < block_count; ++block) {
                        if (!scratch.exponents[block].empty()) {
                            for (std::size_t lane = 0; lane < active_lanes; ++lane) {
                                output[block][first + lane] +=
                                    kernel::template recombine_products<Value>(products[block][lane]);
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
    }    // namespace detail

    /**
     * Reusable baby-step precomputation for Brent-Kung composition with a fixed inner polynomial and divisor. For a
     * block size k, the object stores the reduced powers
     *
     *     1, inner, ..., inner^(k-1) mod B
     *
     * and the giant step inner^k mod B. It can be reused to compose multiple outer polynomials containing at most
     * maximum_outer_coefficient_count coefficients. A divisor context used for later compositions must represent the
     * same polynomial B and contain sufficient inverse precision.
     *
     * Supported Fp12 coefficients may additionally cache their 54 prepared base-field forms. The context's
     * modular_composition_prepared_byte_limit bounds this optional payload; zero, insufficient space, or an
     * unsupported coefficient type keeps the scalar path. Prepared data is owned by this object and read-only
     * during composition. Copies own independent data; scratch belongs to each composition call.
     *
     * @throws std::invalid_argument if maximum_outer_coefficient_count or the configured cached-power limit is zero,
     *         or if reducing inner requires more inverse coefficients than divisor_context stores.
     * @pre inner is a nonempty canonical coefficient polynomial.
     */
    template<detail::SupportsDivrem Backend>
    class polynomial_composition_precomputation {
    public:
        using backend_type = Backend;
        using polynomial_type = typename backend_type::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using coefficient_cache_type = detail::prepared_composition_coefficients<value_type>;

        polynomial_composition_precomputation(
            const polynomial_type &inner, std::size_t maximum_outer_coefficient_count,
            const polynomial_divisor_context<backend_type> &divisor_context,
            polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context) :
            maximum_outer_coefficient_count_(maximum_outer_coefficient_count) {
            if (maximum_outer_coefficient_count_ == 0) {
                throw std::invalid_argument("the maximum outer coefficient count must be positive");
            }

            const std::size_t cached_power_limit = arithmetic_context.options().modular_composition_cached_power_limit;
            if (cached_power_limit == 0) {
                throw std::invalid_argument("the modular-composition cached-power limit must be positive");
            }
            block_size_ = std::min(detail::ceil_sqrt(maximum_outer_coefficient_count_), cached_power_limit);

            if (divisor_context.degree() == 0) {
                return;
            }

            polynomial_type reduced_inner;
            if (inner.size() >= divisor_context.divisor().size()) {
                remainder(reduced_inner, inner, divisor_context, arithmetic_context);
            } else {
                reduced_inner = inner;
            }

            baby_steps_.reserve(block_size_);
            baby_steps_.emplace_back(polynomial_type {value_type::one()});
            if (block_size_ > 1) {
                baby_steps_.emplace_back(std::move(reduced_inner));
            }

            // Build each power from two already cached powers with approximately equal exponents. Even exponents use
            // the dedicated squaring path; odd exponents multiply the floor and ceiling half-powers.
            for (std::size_t exponent = 2; exponent < block_size_; ++exponent) {
                polynomial_type power;
                if (exponent % 2 == 0) {
                    squaremod(power, baby_steps_[exponent / 2], divisor_context, arithmetic_context);
                } else {
                    mulmod(power, baby_steps_[exponent / 2], baby_steps_[exponent / 2 + 1], divisor_context,
                           arithmetic_context);
                }
                baby_steps_.emplace_back(std::move(power));
            }

            if (block_size_ == 1) {
                giant_step_ = std::move(reduced_inner);
            } else if (block_size_ % 2 == 0) {
                squaremod(giant_step_, baby_steps_[block_size_ / 2], divisor_context, arithmetic_context);
            } else {
                mulmod(giant_step_, baby_steps_[block_size_ / 2], baby_steps_[block_size_ / 2 + 1], divisor_context,
                       arithmetic_context);
            }
            if constexpr (coefficient_cache_type::supported) {
                prepared_coefficients_ = coefficient_cache_type(
                    baby_steps_, arithmetic_context.options().modular_composition_prepared_byte_limit,
                    arithmetic_context.options().modular_composition_simd_lanes);
            }
        }

        std::size_t maximum_outer_coefficient_count() const {
            return maximum_outer_coefficient_count_;
        }

        std::size_t block_size() const {
            return block_size_;
        }

        // The composition loop only requests exponents below block_size().
        const polynomial_type &baby_step(std::size_t exponent) const {
            return baby_steps_[exponent];
        }

        const polynomial_type &giant_step() const {
            return giant_step_;
        }

        std::size_t prepared_coefficient_bytes() const {
            if constexpr (coefficient_cache_type::supported) {
                return prepared_coefficients_.bytes();
            } else {
                return 0;
            }
        }

        const coefficient_cache_type &prepared_coefficients() const {
            return prepared_coefficients_;
        }

    private:
        std::size_t maximum_outer_coefficient_count_;
        std::size_t block_size_ = 0;
        std::vector<polynomial_type> baby_steps_;
        polynomial_type giant_step_;
        [[no_unique_address]] coefficient_cache_type prepared_coefficients_;
    };

    namespace detail {
        // Add one evaluated outer block to an existing accumulator. The cache must have a nonconstant divisor;
        // outer must fit it, block_index must exist, and result must not alias outer or the cached powers.
        template<SupportsDivrem Backend>
        void add_composition_block(typename Backend::polynomial_type &result,
                                   const typename Backend::polynomial_type &outer,
                                   const polynomial_composition_precomputation<Backend> &precomputation,
                                   std::size_t block_index) {
            using polynomial_type = typename Backend::polynomial_type;
            using value_type = typename polynomial_type::value_type;
            using cache_type = typename polynomial_composition_precomputation<Backend>::coefficient_cache_type;
            if constexpr (cache_type::supported) {
                if (precomputation.prepared_coefficients().use_for_block(outer, block_index)) {
                    typename cache_type::workspace scratch;
                    std::array<polynomial_type, cache_type::block_tile_size> blocks;
                    precomputation.prepared_coefficients().evaluate_blocks(blocks, outer, block_index, 1, scratch);
                    addition(result, result, blocks[0]);
                    return;
                }
            }

            const std::size_t block_size = precomputation.block_size();
            const std::size_t first = block_index * block_size;
            const std::size_t count = std::min(block_size, outer.size() - first);
            for (std::size_t exponent = 0; exponent < count; ++exponent) {
                const value_type &scalar = outer[first + exponent];
                if (scalar == value_type {}) {
                    continue;
                }
                if (exponent == 0) {
                    result[0] = result[0] + scalar;
                    continue;
                }
                const polynomial_type &power = precomputation.baby_step(exponent);
                if (result.size() < power.size()) {
                    result.resize(power.size(), value_type {});
                }
                for (std::size_t coefficient = 0; coefficient < power.size(); ++coefficient) {
                    result[coefficient] = result[coefficient] + power[coefficient] * scalar;
                }
            }
            condense(result);
        }
    }    // namespace detail

    /**
     * Compute output = outer(inner(X)) mod B by Horner's rule, where B is the nonzero polynomial stored in
     * divisor_context. This Horner algorithm performs one modular multiplication for each coefficient of outer except
     * its leading coefficient and is intended as a correctness oracle for faster modular composition algorithms. The
     * result is canonical, has degree less than degree(B) when B is nonconstant, and may alias either input.
     *
     * The inner polynomial is reduced before composition, and every Horner product is reduced before the next
     * coefficient is added. For d = degree(B), the products of reduced polynomials require at most d - 1 inverse
     * coefficients. Reducing an unreduced inner polynomial may require greater precomputed precision. Modulo a
     * nonzero constant, every composition is the zero polynomial.
     *
     * @throws std::invalid_argument if the precomputed inverse has insufficient precision.
     * @pre outer and inner are nonempty canonical coefficient polynomials.
     */
    template<detail::SupportsDivrem Backend>
    void compose_mod_reference(typename Backend::polynomial_type &output,
                               const typename Backend::polynomial_type &outer,
                               const typename Backend::polynomial_type &inner,
                               const polynomial_divisor_context<Backend> &divisor_context,
                               polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if (divisor_context.degree() == 0) {
            output.assign(1, value_type {});
            return;
        }
        if (outer.size() == 1) {
            output = outer;
            return;
        }

        // Reference an already reduced inner polynomial directly to avoid copying it; use local storage only when
        // reduction is required.
        polynomial_type reduced_inner_storage;
        const polynomial_type *reduced_inner = &inner;
        if (inner.size() >= divisor_context.divisor().size()) {
            remainder(reduced_inner_storage, inner, divisor_context, arithmetic_context);
            reduced_inner = &reduced_inner_storage;
        }

        polynomial_type result {outer.back()};
        for (std::size_t i = outer.size() - 1; i-- > 0;) {
            mulmod(result, result, *reduced_inner, divisor_context, arithmetic_context);
            result[0] = result[0] + outer[i];
        }

        output = std::move(result);
    }

    /**
     * Compute output = outer(inner(X)) mod B using a reusable Brent-Kung precomputation for inner and B. The result is
     * canonical, has degree less than degree(B) when B is nonconstant, and may alias outer.
     *
     * For a block size k, split outer into polynomials F_j of degree less than k:
     *
     *     outer(Y) = F_0(Y) + F_1(Y) * Y^k + F_2(Y) * Y^(2k) + ...
     *
     * The precomputation caches the baby steps 1, inner, ..., inner^(k-1) modulo B. Each F_j(inner) is then a linear
     * combination of these cached powers. It also stores the giant step inner^k modulo B, used to combine the block
     * values by Horner's rule:
     *
     *     (...(F_last(inner) * giant_step + F_previous(inner)) * giant_step + ...) + F_0(inner).
     *
     * Constructing the precomputation uses about k modular multiplications, and each composition uses about L/k more
     * when outer has L coefficients. Thus a one-off composition uses about k + L/k modular multiplications, minimized
     * by k = ceil(sqrt(L)). Forming the block values also requires O(L * degree(B)) coefficient operations; this
     * quadratic work can dominate when L and degree(B) are both large. The cached-power limit may select a smaller k,
     * storing at most k * degree(B) coefficients at the cost of more giant-step multiplications.
     *
     * When the precomputation contains prepared coefficient forms, up to four adjacent blocks with sufficiently
     * many nonzero nonconstant terms are evaluated together using fused base-field dot products. Short or sparse
     * blocks use scalar multiplication. Both paths preserve the descending Horner order. At most four block
     * polynomials are buffered; each selected outer coefficient is prepared once, inside the block kernel.
     *
     * @throws std::invalid_argument if outer exceeds the maximum coefficient count selected when the precomputation
     *         was constructed.
     * @pre outer is a nonempty canonical coefficient polynomial. divisor_context represents the same divisor used to
     *      construct precomputation and contains sufficient inverse precision.
     */
    template<detail::SupportsDivrem Backend>
    void compose_mod(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &outer,
                     const polynomial_composition_precomputation<Backend> &precomputation,
                     const polynomial_divisor_context<Backend> &divisor_context,
                     polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;
        using value_type = typename polynomial_type::value_type;

        if (outer.size() > precomputation.maximum_outer_coefficient_count()) {
            throw std::invalid_argument("outer polynomial exceeds the modular-composition precomputation");
        }
        if (divisor_context.degree() == 0) {
            output.assign(1, value_type {});
            return;
        }
        if (outer.size() == 1) {
            output = outer;
            return;
        }

        const std::size_t block_size = precomputation.block_size();
        const std::size_t block_count = 1 + (outer.size() - 1) / block_size;
        using cache_type = typename polynomial_composition_precomputation<Backend>::coefficient_cache_type;
        if constexpr (cache_type::supported) {
            if (precomputation.prepared_coefficients().available()) {
                typename cache_type::workspace scratch;
                std::array<polynomial_type, cache_type::block_tile_size> blocks;
                polynomial_type result;
                bool leading = true;
                for (std::size_t end = block_count; end > 0;) {
                    std::size_t first = end;
                    while (first > 0 && end - first < cache_type::block_tile_size &&
                           precomputation.prepared_coefficients().use_for_block(outer, first - 1)) {
                        --first;
                    }
                    if (first == end) {
                        if (!leading) {
                            mulmod(result, result, precomputation.giant_step(), divisor_context, arithmetic_context);
                        }
                        detail::add_composition_block(result, outer, precomputation, end - 1);
                        leading = false;
                        --end;
                        continue;
                    }
                    const std::size_t count = end - first;
                    precomputation.prepared_coefficients().evaluate_blocks(blocks, outer, first, count, scratch);
                    // Evaluate blocks together, then consume them in the original descending Horner order.
                    for (std::size_t block = count; block-- > 0;) {
                        if (leading) {
                            result = std::move(blocks[block]);
                            leading = false;
                        } else {
                            mulmod(result, result, precomputation.giant_step(), divisor_context, arithmetic_context);
                            addition(result, result, blocks[block]);
                        }
                    }
                    end = first;
                }
                output = std::move(result);
                return;
            }
        }
        polynomial_type result;
        detail::add_composition_block(result, outer, precomputation, block_count - 1);
        for (std::size_t block = block_count - 1; block-- > 0;) {
            mulmod(result, result, precomputation.giant_step(), divisor_context, arithmetic_context);
            detail::add_composition_block(result, outer, precomputation, block);
        }

        output = std::move(result);
    }

    /**
     * Compute output = outer(inner(X)) mod B using blocked Brent-Kung composition, where B is the nonzero polynomial
     * stored in divisor_context. This overload constructs the baby-step precomputation for one composition; callers
     * composing several outer polynomials with the same inner and divisor should construct and reuse a
     * polynomial_composition_precomputation instead.
     *
     * The result is canonical, has degree less than degree(B) when B is nonconstant, and may alias either input. The
     * inner polynomial is reduced before its powers are constructed. Reducing an unreduced inner polynomial may
     * require more precomputed inverse coefficients than products of reduced operands.
     *
     * @throws std::invalid_argument if the cached-power limit is zero or the precomputed inverse has insufficient
     *         precision.
     * @pre outer and inner are nonempty canonical coefficient polynomials.
     */
    template<detail::SupportsDivrem Backend>
    void compose_mod(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &outer,
                     const typename Backend::polynomial_type &inner,
                     const polynomial_divisor_context<Backend> &divisor_context,
                     polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using value_type = typename Backend::polynomial_type::value_type;

        if (arithmetic_context.options().modular_composition_cached_power_limit == 0) {
            throw std::invalid_argument("the modular-composition cached-power limit must be positive");
        }
        if (divisor_context.degree() == 0) {
            output.assign(1, value_type {});
            return;
        }
        if (outer.size() == 1) {
            output = outer;
            return;
        }

        polynomial_composition_precomputation<Backend> precomputation(inner, outer.size(), divisor_context,
                                                                      arithmetic_context);
        compose_mod(output, outer, precomputation, divisor_context, arithmetic_context);
    }

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_COMPOSITION_HPP
