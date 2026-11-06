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
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/math/detail/integer_sqrt.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/detail/prepared_composition_coefficients.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_modular_arithmetic.hpp>

namespace nil::crypto3::math {

    namespace detail {
        template<typename Backend>
        using product_accumulator_preparation =
            decltype(std::declval<polynomial_arithmetic::polynomial_context<Backend> &>().make_product_accumulator(
                std::size_t {}, std::size_t {}));

        template<typename Backend>
        concept SupportsPreparedProductSum =
            SupportsPreparedLowProduct<Backend> &&
            requires(polynomial_arithmetic::polynomial_context<Backend> &context,
                     typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                     typename product_accumulator_preparation<Backend>::value_type &accumulator, std::size_t count) {
                { context.make_product_accumulator(count, count).has_value() } -> std::same_as<bool>;
                {
                    context.make_product_accumulator(count, count)->required_spectrum_bytes()
                } -> std::convertible_to<std::size_t>;
                {
                    context.prepare_low_product(input, count, count)->storage_bytes()
                } -> std::convertible_to<std::size_t>;
                {
                    context.try_accumulate_prepared_product(accumulator, input,
                                                            *context.prepare_low_product(input, count, count))
                } -> std::same_as<bool>;
                { context.try_finish_product_sum(output, accumulator) } -> std::same_as<bool>;
            };

        // Follow divisor_prepared_products: backends without the optional capability carry no state.
        template<typename Backend>
        struct composition_prepared_giant_steps { };

        template<SupportsPreparedProductSum Backend>
        struct composition_prepared_giant_steps<Backend> {
            using context_type = polynomial_arithmetic::polynomial_context<Backend>;
            using preparation_type = decltype(std::declval<context_type &>().prepare_low_product(
                std::declval<const typename Backend::polynomial_type &>(), std::size_t {}, std::size_t {}));
            using accumulator_type = product_accumulator_preparation<Backend>;
            std::vector<typename preparation_type::value_type> powers;
            // This seed stays empty. Each composition copies its certified CRT setup into private scratch.
            accumulator_type accumulator;
            std::size_t bytes = 0;
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
     * The default block size is ceil(sqrt(maximum_outer_coefficient_count)). An explicit context option may select
     * a larger block for repeated compositions, subject to the cached-power limit and maximum outer size.
     *
     * Supported Fp12 coefficients may additionally cache prepared base-field forms (23 for BN254 SIMD,
     * 54 for the scalar tower kernel). The context's
     * modular_composition_prepared_byte_limit bounds this optional payload; zero, insufficient space, or an
     * unsupported coefficient type keeps the scalar path. Prepared data is owned by this object and read-only
     * during composition. Copies own independent data; scratch belongs to each composition call.
     * An optional constructor byte limit further caps the prepared coefficient payload without changing the
     * arithmetic context or other caches. It does not bound ordinary cached powers or per-call scratch.
     *
     * A separately enabled giant-step group cache prepares G, G^2, ..., G^g modulo B, where G=inner^k mod B.
     * modular_composition_giant_step_prepared_byte_limit bounds those spectra independently of the coefficient
     * cache. Group size zero, insufficient space/capacity, or an unsupported backend keeps Horner evaluation.
     * Group size one reuses G's transform without changing the number of modular reductions.
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
        using giant_cache_type = detail::composition_prepared_giant_steps<backend_type>;

        polynomial_composition_precomputation(
            const polynomial_type &inner, std::size_t maximum_outer_coefficient_count,
            const polynomial_divisor_context<backend_type> &divisor_context,
            polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context,
            std::size_t prepared_byte_limit = std::numeric_limits<std::size_t>::max()) :
            maximum_outer_coefficient_count_(maximum_outer_coefficient_count) {
            if (maximum_outer_coefficient_count_ == 0) {
                throw std::invalid_argument("the maximum outer coefficient count must be positive");
            }

            const std::size_t cached_power_limit = arithmetic_context.options().modular_composition_cached_power_limit;
            if (cached_power_limit == 0) {
                throw std::invalid_argument("the modular-composition cached-power limit must be positive");
            }
            const std::size_t requested_block_size = arithmetic_context.options().modular_composition_block_size;
            const std::size_t selected_block_size =
                requested_block_size == 0 ? detail::ceil_sqrt(maximum_outer_coefficient_count_) :
                                            std::min(requested_block_size, maximum_outer_coefficient_count_);
            block_size_ = std::min(selected_block_size, cached_power_limit);

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
                    baby_steps_,
                    std::min(prepared_byte_limit, arithmetic_context.options().modular_composition_prepared_byte_limit),
                    arithmetic_context.options().modular_composition_simd_lanes);
            }
            prepare_giant_steps(divisor_context, arithmetic_context);
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

        std::size_t giant_step_group_size() const {
            if constexpr (detail::SupportsPreparedProductSum<backend_type>) {
                return prepared_giant_steps_.powers.size();
            } else {
                return 0;
            }
        }

        std::size_t prepared_giant_step_bytes() const {
            if constexpr (detail::SupportsPreparedProductSum<backend_type>) {
                return prepared_giant_steps_.bytes;
            } else {
                return 0;
            }
        }

        const giant_cache_type &prepared_giant_steps() const {
            return prepared_giant_steps_;
        }

    private:
        void prepare_giant_steps(const polynomial_divisor_context<backend_type> &divisor_context,
                                 polynomial_arithmetic::polynomial_context<backend_type> &arithmetic_context) {
            if constexpr (detail::SupportsPreparedProductSum<backend_type>) {
                const auto &options = arithmetic_context.options();
                const std::size_t group = std::min(options.modular_composition_giant_step_group_size,
                                                   (maximum_outer_coefficient_count_ - 1) / block_size_);
                const std::size_t limit = options.modular_composition_giant_step_prepared_byte_limit;
                const std::size_t degree = divisor_context.degree();
                if (group == 0 || limit == 0 || degree > std::numeric_limits<std::size_t>::max() / 2) {
                    return;
                }
                giant_cache_type prepared;
                prepared.accumulator = arithmetic_context.make_product_accumulator(degree, group);
                if (!prepared.accumulator) {
                    return;
                }
                const std::size_t spectrum_bytes = prepared.accumulator->required_spectrum_bytes();
                if (spectrum_bytes == 0 || group > limit / spectrum_bytes) {
                    return;
                }
                prepared.powers.reserve(group);
                polynomial_type power = giant_step_;
                for (std::size_t exponent = 1; exponent <= group; ++exponent) {
                    if (exponent > 1) {
                        mulmod(power, power, giant_step_, divisor_context, arithmetic_context);
                    }
                    // Every fixed spectrum uses the group's full transform, including sparse or constant powers.
                    auto padded = power;
                    padded.resize(degree, value_type::zero());
                    auto operand = arithmetic_context.prepare_low_product(padded, degree, 2 * degree - 1);
                    if (!operand || operand->storage_bytes() > limit - prepared.bytes) {
                        return;
                    }
                    prepared.bytes += operand->storage_bytes();
                    prepared.powers.push_back(std::move(*operand));
                }
                prepared_giant_steps_ = std::move(prepared);
            }
        }

        std::size_t maximum_outer_coefficient_count_;
        std::size_t block_size_ = 0;
        std::vector<polynomial_type> baby_steps_;
        polynomial_type giant_step_;
        [[no_unique_address]] coefficient_cache_type prepared_coefficients_;
        [[no_unique_address]] giant_cache_type prepared_giant_steps_;
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

        // Consume the same scalar/prepared block values as Horner, but combine g steps at a time:
        // A*G^g + C_(g-1)*G^(g-1) + ... + C_1*G + C_0, with every cached power reduced modulo B.
        // The sum has degree <= 2*degree(B)-2 regardless of g, so one reduction suffices per group.
        // Only publish the final output: a declined backend operation can restart ordinary Horner safely.
        template<SupportsDivrem Backend>
            requires SupportsPreparedProductSum<Backend>
        bool try_compose_mod_grouped(typename Backend::polynomial_type &output,
                                     const typename Backend::polynomial_type &outer,
                                     const polynomial_composition_precomputation<Backend> &precomputation,
                                     const polynomial_divisor_context<Backend> &divisor_context,
                                     polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
            using polynomial_type = typename Backend::polynomial_type;
            using cache_type = typename polynomial_composition_precomputation<Backend>::coefficient_cache_type;
            const auto &giants = precomputation.prepared_giant_steps();
            const std::size_t group = giants.powers.size();
            if (group == 0 || !giants.accumulator) {
                return false;
            }
            auto accumulator = *giants.accumulator;
            polynomial_type result, product;
            bool has_products = false;
            const auto append = [&](const polynomial_type &input, std::size_t exponent) {
                if (!arithmetic_context.try_accumulate_prepared_product(accumulator, input,
                                                                        giants.powers[exponent - 1])) {
                    return false;
                }
                has_products = true;
                return true;
            };
            const auto consume = [&](polynomial_type &block, std::size_t index) {
                const std::size_t exponent = index % group;
                if (exponent != 0) {
                    return append(block, exponent);
                }
                if (has_products) {
                    if (!arithmetic_context.try_finish_product_sum(product, accumulator)) {
                        return false;
                    }
                    remainder(result, product, divisor_context, arithmetic_context);
                    addition(result, result, block);
                } else {
                    result = std::move(block);
                }
                has_products = false;
                return index == 0 || append(result, group);
            };
            const auto consume_scalar = [&](std::size_t index) {
                polynomial_type block;
                add_composition_block(block, outer, precomputation, index);
                return consume(block, index);
            };
            std::size_t end = 1 + (outer.size() - 1) / precomputation.block_size();
            if constexpr (cache_type::supported) {
                if (precomputation.prepared_coefficients().available()) {
                    typename cache_type::workspace scratch;
                    std::array<polynomial_type, cache_type::block_tile_size> blocks;
                    while (end > 0) {
                        std::size_t first = end;
                        while (first > 0 && end - first < cache_type::block_tile_size &&
                               precomputation.prepared_coefficients().use_for_block(outer, first - 1)) {
                            --first;
                        }
                        if (first == end) {
                            if (!consume_scalar(--end)) {
                                return false;
                            }
                            continue;
                        }
                        const std::size_t count = end - first;
                        precomputation.prepared_coefficients().evaluate_blocks(blocks, outer, first, count, scratch);
                        for (std::size_t block = count; block-- > 0;) {
                            if (!consume(blocks[block], first + block)) {
                                return false;
                            }
                        }
                        end = first;
                    }
                }
            }
            while (end > 0) {
                if (!consume_scalar(--end)) {
                    return false;
                }
            }
            output = std::move(result);
            return true;
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
     * blocks use scalar multiplication. Both paths visit blocks in descending order. At most four block
     * polynomials are buffered; each selected outer coefficient is prepared once, inside the block kernel.
     * An optional prepared giant-step group combines several products in transform space, then reconstructs and
     * reduces their sum once. Unsupported or declined grouped operations fall back to the ordinary Horner loop.
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
        if constexpr (detail::SupportsPreparedProductSum<Backend>) {
            if (block_count > 1 && precomputation.giant_step_group_size() != 0 &&
                detail::try_compose_mod_grouped(output, outer, precomputation, divisor_context, arithmetic_context)) {
                return;
            }
        }
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
