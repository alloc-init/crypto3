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

#define BOOST_TEST_MODULE polynomial_composition_test

#include <cstddef>
#include <limits>
#include <optional>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_composition.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace polynomial_arithmetic = math::polynomial_arithmetic;
    namespace fields = nil::crypto3::algebra::fields;

    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq_value_type = fq_field_type::value_type;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using fq12_value_type = fq12_field_type::value_type;

    fq12_value_type fq12_value(std::size_t first_coordinate) {
        fq12_value_type value = fq12_value_type::zero();
        for (std::size_t i = 0; i < fq12_field_type::arity; ++i) {
            value.coordinate(i) = fq_value_type(first_coordinate + i);
        }
        return value;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(polynomial_composition_test_suite)

BOOST_AUTO_TEST_CASE(reference_composition_matches_a_known_result_and_may_alias_inputs) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    // Modulo X^2 + 1, (1 + 2Y + 3Y^2) evaluated at Y = X is -2 + 2X.
    const polynomial_type outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3)};
    const polynomial_type inner = {fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    const polynomial_type expected = {-fq_value_type(2), fq_value_type(2)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);

    polynomial_type result;
    math::compose_mod_reference(result, outer, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == expected);

    polynomial_type outer_alias = outer;
    math::compose_mod_reference(outer_alias, outer_alias, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(outer_alias == expected);

    polynomial_type inner_alias = inner;
    math::compose_mod_reference(inner_alias, outer, inner_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(inner_alias == expected);
}

BOOST_AUTO_TEST_CASE(reference_composition_reduces_the_inner_polynomial) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type outer = {fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type inner = {fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8)};
    const polynomial_type divisor = {fq_value_type(1), fq_value_type(1), fq_value_type(1)};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 2, arithmetic_context);

    polynomial_type reduced_inner;
    math::remainder(reduced_inner, inner, divisor_context, arithmetic_context);

    polynomial_type from_unreduced;
    polynomial_type from_reduced;
    math::compose_mod_reference(from_unreduced, outer, inner, divisor_context, arithmetic_context);
    math::compose_mod_reference(from_reduced, outer, reduced_inner, divisor_context, arithmetic_context);
    BOOST_CHECK(from_unreduced == from_reduced);

    polynomial_type inner_alias = inner;
    math::compose_mod_reference(inner_alias, outer, inner_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(inner_alias == from_unreduced);
}

BOOST_AUTO_TEST_CASE(composition_handles_zero_constant_and_constant_modulus_cases) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    const polynomial_type divisor = {fq_value_type(2), fq_value_type(3), fq_value_type(1)};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);
    const polynomial_type inner = {fq_value_type(4), fq_value_type(5), fq_value_type(6), fq_value_type(7)};

    polynomial_type result;
    math::compose_mod_reference(result, polynomial_type {fq_value_type::zero()}, inner, divisor_context,
                                arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
    math::compose_mod(result, polynomial_type {fq_value_type::zero()}, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));

    math::compose_mod_reference(result, polynomial_type {fq_value_type(7)}, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type(7)}));
    math::compose_mod(result, polynomial_type {fq_value_type(7)}, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type(7)}));

    const polynomial_type constant_divisor = {fq_value_type(11)};
    math::polynomial_divisor_context<backend_type> constant_context(constant_divisor, 1, arithmetic_context);
    math::compose_mod_reference(result, polynomial_type {fq_value_type(1), fq_value_type(2)}, inner, constant_context,
                                arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
    math::compose_mod(result, polynomial_type {fq_value_type(1), fq_value_type(2)}, inner, constant_context,
                      arithmetic_context);
    BOOST_CHECK(result == polynomial_type({fq_value_type::zero()}));
}

BOOST_AUTO_TEST_CASE(composition_rejects_insufficient_inverse_precision_without_changing_output) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);

    const polynomial_type divisor = {fq_value_type::one(), fq_value_type::zero(), fq_value_type::one()};
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 1, arithmetic_context);
    const polynomial_type outer = {fq_value_type(2), fq_value_type(3)};
    const polynomial_type inner = {fq_value_type(4), fq_value_type(5), fq_value_type(6), fq_value_type(7)};

    polynomial_type output = {fq_value_type(17)};
    BOOST_CHECK_THROW(math::compose_mod_reference(output, outer, inner, divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(output == polynomial_type({fq_value_type(17)}));

    BOOST_CHECK_THROW(math::compose_mod(output, outer, inner, divisor_context, arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(output == polynomial_type({fq_value_type(17)}));
}

BOOST_AUTO_TEST_CASE(brent_kung_composition_matches_the_reference_at_different_memory_limits) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3), fq_value_type(4),
                                   fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8),
                                   fq_value_type(9), fq_value_type(10)};
    const polynomial_type inner = {fq_value_type(11), fq_value_type(12), fq_value_type(13),
                                   fq_value_type(14), fq_value_type(15), fq_value_type(16)};
    const polynomial_type divisor = {fq_value_type(17), fq_value_type(18), fq_value_type(19),   fq_value_type(20),
                                     fq_value_type(21), fq_value_type(22), fq_value_type::one()};

    polynomial_arithmetic::polynomial_context_options reference_options;
    reference_options.basecase_divisor_coefficient_cutoff = 0;
    reference_options.basecase_quotient_coefficient_cutoff = 0;
    polynomial_arithmetic::polynomial_context<backend_type> reference_context(backend_type {}, reference_options);
    math::polynomial_divisor_context<backend_type> reference_divisor(divisor, 5, reference_context);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, reference_divisor, reference_context);

    for (const std::size_t cached_power_limit : {1, 2, 3, 4, 100}) {
        polynomial_arithmetic::polynomial_context_options options = reference_options;
        options.modular_composition_cached_power_limit = cached_power_limit;
        polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> divisor_context(divisor, 5, arithmetic_context);

        polynomial_type result;
        math::compose_mod(result, outer, inner, divisor_context, arithmetic_context);
        BOOST_CHECK(result == expected);
    }
}

BOOST_AUTO_TEST_CASE(brent_kung_precomputation_is_reusable_and_enforces_its_outer_size) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type first_outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3), fq_value_type(4),
                                         fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8),
                                         fq_value_type(9), fq_value_type(10)};
    const polynomial_type second_outer = {fq_value_type(11), fq_value_type(12), fq_value_type(13),
                                          fq_value_type(14), fq_value_type(15), fq_value_type(16)};
    const polynomial_type inner = {fq_value_type(17), fq_value_type(18), fq_value_type(19),
                                   fq_value_type(20), fq_value_type(21), fq_value_type(22)};
    const polynomial_type divisor = {fq_value_type(23), fq_value_type(24), fq_value_type(25),   fq_value_type(26),
                                     fq_value_type(27), fq_value_type(28), fq_value_type::one()};

    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context;
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 5, arithmetic_context);
    math::polynomial_composition_precomputation<backend_type> precomputation(inner, first_outer.size(), divisor_context,
                                                                             arithmetic_context);

    for (const polynomial_type &outer : {first_outer, second_outer}) {
        polynomial_type expected;
        math::compose_mod_reference(expected, outer, inner, divisor_context, arithmetic_context);

        polynomial_type result;
        math::compose_mod(result, outer, precomputation, divisor_context, arithmetic_context);
        BOOST_CHECK(result == expected);

        polynomial_type aliased = outer;
        math::compose_mod(aliased, aliased, precomputation, divisor_context, arithmetic_context);
        BOOST_CHECK(aliased == expected);
    }

    polynomial_type oversized_outer = first_outer;
    oversized_outer.push_back(fq_value_type(11));
    polynomial_type unchanged = {fq_value_type(17)};
    BOOST_CHECK_THROW(
        math::compose_mod(unchanged, oversized_outer, precomputation, divisor_context, arithmetic_context),
        std::invalid_argument);
    BOOST_CHECK(unchanged == polynomial_type({fq_value_type(17)}));

    const polynomial_type constant_divisor = {fq_value_type(29)};
    math::polynomial_divisor_context<backend_type> constant_divisor_context(constant_divisor, 1, arithmetic_context);
    math::polynomial_composition_precomputation<backend_type> constant_precomputation(
        inner, first_outer.size(), constant_divisor_context, arithmetic_context);
    BOOST_CHECK_THROW(math::compose_mod(unchanged, oversized_outer, constant_precomputation, constant_divisor_context,
                                        arithmetic_context),
                      std::invalid_argument);
    BOOST_CHECK(unchanged == polynomial_type({fq_value_type(17)}));

    BOOST_CHECK_THROW(
        math::polynomial_composition_precomputation<backend_type>(inner, 0, divisor_context, arithmetic_context),
        std::invalid_argument);

    polynomial_arithmetic::polynomial_context_options invalid_options;
    invalid_options.modular_composition_cached_power_limit = 0;
    polynomial_arithmetic::polynomial_context<backend_type> invalid_arithmetic_context(backend_type {},
                                                                                       invalid_options);
    math::polynomial_divisor_context<backend_type> invalid_divisor_context(divisor, 5, invalid_arithmetic_context);
    BOOST_CHECK_THROW(math::polynomial_composition_precomputation<backend_type>(
                          inner, first_outer.size(), invalid_divisor_context, invalid_arithmetic_context),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(brent_kung_composition_may_alias_inputs_and_rejects_a_zero_cache_limit) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;

    const polynomial_type outer = {fq_value_type(1), fq_value_type(2), fq_value_type(3), fq_value_type(4)};
    const polynomial_type inner = {fq_value_type(5), fq_value_type(6), fq_value_type(7), fq_value_type(8)};
    const polynomial_type divisor = {fq_value_type(9), fq_value_type(10), fq_value_type(11), fq_value_type::one()};

    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_cached_power_limit = 2;
    polynomial_arithmetic::polynomial_context<backend_type> arithmetic_context(backend_type {}, options);
    math::polynomial_divisor_context<backend_type> divisor_context(divisor, 2, arithmetic_context);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, divisor_context, arithmetic_context);

    polynomial_type outer_alias = outer;
    math::compose_mod(outer_alias, outer_alias, inner, divisor_context, arithmetic_context);
    BOOST_CHECK(outer_alias == expected);

    polynomial_type inner_alias = inner;
    math::compose_mod(inner_alias, outer, inner_alias, divisor_context, arithmetic_context);
    BOOST_CHECK(inner_alias == expected);

    options.modular_composition_cached_power_limit = 0;
    polynomial_arithmetic::polynomial_context<backend_type> invalid_context(backend_type {}, options);
    math::polynomial_divisor_context<backend_type> invalid_divisor(divisor, 2, invalid_context);
    polynomial_type output = {fq_value_type(17)};
    BOOST_CHECK_THROW(math::compose_mod(output, outer, inner, invalid_divisor, invalid_context), std::invalid_argument);
    BOOST_CHECK(output == polynomial_type({fq_value_type(17)}));
}

BOOST_AUTO_TEST_CASE(composition_matches_between_schoolbook_and_mixed_radix_backends) {
    using schoolbook_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using mixed_radix_backend = polynomial_arithmetic::mixed_radix_backend<fq_field_type, fq12_value_type>;
    using polynomial_type = typename schoolbook_backend::polynomial_type;

    const polynomial_type outer = {fq12_value(1), fq12_value(13), fq12_value(25), fq12_value(37)};
    const polynomial_type inner = {fq12_value(49), fq12_value(61), fq12_value(73)};
    const polynomial_type divisor = {fq12_value(85), fq12_value(97), fq12_value(109), fq12_value(121)};

    polynomial_arithmetic::polynomial_context<schoolbook_backend> schoolbook_context;
    math::polynomial_divisor_context<schoolbook_backend> schoolbook_divisor(divisor, 2, schoolbook_context);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, schoolbook_divisor, schoolbook_context);

    polynomial_arithmetic::polynomial_context_options options;
    options.basecase_divisor_coefficient_cutoff = 0;
    options.basecase_quotient_coefficient_cutoff = 0;
    options.modular_composition_cached_power_limit = 2;
    polynomial_arithmetic::polynomial_context<mixed_radix_backend> mixed_radix_context(mixed_radix_backend(18),
                                                                                       options);
    math::polynomial_divisor_context<mixed_radix_backend> mixed_radix_divisor(divisor, 2, mixed_radix_context);
    polynomial_type result;
    math::compose_mod_reference(result, outer, inner, mixed_radix_divisor, mixed_radix_context);
    BOOST_CHECK(result == expected);

    math::compose_mod(result, outer, inner, mixed_radix_divisor, mixed_radix_context);
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_CASE(prepared_composition_cache_obeys_its_budget_and_pads_short_powers) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using precomputation_type = math::polynomial_composition_precomputation<backend_type>;
    using kernel = fields::detail::fp12_fast::fp12_fast<fq12_field_type::extension_policy::fast_params>;
    // k=5, inner=X: the nonconstant baby steps have lengths 2,3,4,5. Padding is included in the budget.
    const std::size_t required = 4 * 5 * sizeof(kernel::prepared_product_operand);
    const polynomial_type inner = {fq12_value_type::zero(), fq12_value_type::one()};
    polynomial_type divisor(8, fq12_value_type::zero());
    divisor[0] = fq12_value(3);
    divisor.back() = fq12_value(17);    // Nonmonic divisor.
    polynomial_type outer(21);
    for (std::size_t i = 0; i < outer.size(); ++i) {
        outer[i] = fq12_value(23 + 12 * i);
    }
    for (std::size_t limit :
         {std::size_t(0), required - 1, required, required + 1, std::numeric_limits<std::size_t>::max()}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = 5;
        options.modular_composition_prepared_byte_limit = limit;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, 6, context);
        const precomputation_type cache(inner, outer.size(), division, context);
        BOOST_CHECK_EQUAL(cache.block_size(), 5);
        BOOST_CHECK_EQUAL(cache.prepared_coefficient_bytes(), limit >= required ? required : 0);
        polynomial_type expected, actual;
        math::compose_mod_reference(expected, outer, inner, division, context);
        math::compose_mod(actual, outer, cache, division, context);
        BOOST_CHECK(actual == expected);

        // Check the shared single-block entry point too, including an existing accumulator and cancellation.
        for (std::size_t block = 0; block < 5; ++block) {
            const std::size_t first = block * 5;
            const std::size_t end = std::min(first + 5, outer.size());
            polynomial_type slice(outer.begin() + first, outer.begin() + end), evaluated;
            math::compose_mod_reference(evaluated, slice, inner, division, context);
            actual = {fq12_value(77), fq12_value(91)};
            math::addition(expected, actual, evaluated);
            math::detail::add_composition_block(actual, outer, cache, block);
            BOOST_CHECK(actual == expected);
            math::subtraction(actual, polynomial_type {}, evaluated);
            math::detail::add_composition_block(actual, outer, cache, block);
            BOOST_CHECK(actual == polynomial_type {fq12_value_type::zero()});
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_cache_owns_its_data_and_survives_copy_and_move) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    using precomputation_type = math::polynomial_composition_precomputation<backend_type>;
    const polynomial_type inner = {fq12_value(1), fq12_value(13), fq12_value(25)};
    const polynomial_type divisor = {fq12_value(37), fq12_value(49), fq12_value(61), fq12_value(73), fq12_value(85)};
    std::optional<precomputation_type> saved;
    {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_prepared_byte_limit = 1024 * 1024;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, 3, context);
        polynomial_type temporary_inner = inner;
        saved.emplace(temporary_inner, 25, division, context);
        temporary_inner.assign(1, fq12_value(999));
    }
    BOOST_REQUIRE_GT(saved->prepared_coefficient_bytes(), 0);
    auto copy = *saved;
    const std::size_t bytes = copy.prepared_coefficient_bytes();
    saved.reset();
    auto moved = std::move(copy);
    BOOST_CHECK_EQUAL(moved.prepared_coefficient_bytes(), bytes);
    polynomial_arithmetic::polynomial_context<backend_type> context;
    math::polynomial_divisor_context<backend_type> division(divisor, 3, context);
    for (std::size_t count : {1, 4, 5, 6, 19, 20, 21, 25}) {
        polynomial_type outer(count);
        for (std::size_t i = 0; i < count; ++i) {
            outer[i] = i % 3 == 0 ? fq12_value_type::zero() : fq12_value(i * 12 + count);
        }
        outer.back() = fq12_value(123);
        polynomial_type expected, actual;
        math::compose_mod_reference(expected, outer, inner, division, context);
        math::compose_mod(actual, outer, moved, division, context);
        BOOST_CHECK(actual == expected);
        math::compose_mod(outer, outer, moved, division, context);
        BOOST_CHECK(outer == expected);
    }
    polynomial_type oversized(26, fq12_value(1)), unchanged {fq12_value(91)};
    BOOST_CHECK_THROW(math::compose_mod(unchanged, oversized, moved, division, context), std::invalid_argument);
    BOOST_CHECK(unchanged == polynomial_type {fq12_value(91)});
}

BOOST_AUTO_TEST_CASE(prepared_composition_dispatch_counts_nonzero_nonconstant_terms) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    constexpr std::size_t block_size = 17;
    polynomial_type divisor(10), inner(9);
    for (std::size_t i = 0; i < divisor.size(); ++i)
        divisor[i] = fq12_value(12 * i + 1);
    for (std::size_t i = 0; i < inner.size(); ++i)
        inner[i] = fq12_value(12 * i + 121);

    for (std::size_t lanes : {0, 4, 8}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = block_size;
        options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
        options.modular_composition_simd_lanes = lanes;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, 8, context);
        const math::polynomial_composition_precomputation<backend_type> cache(inner, block_size * block_size, division,
                                                                              context);
        BOOST_REQUIRE_EQUAL(cache.block_size(), block_size);
        BOOST_REQUIRE(cache.prepared_coefficients().available());
#if defined(CRYPTO3_HAS_AVX512IFMA)
        BOOST_CHECK_EQUAL(cache.prepared_coefficients().simd_lanes(), lanes);
#else
        BOOST_CHECK_EQUAL(cache.prepared_coefficients().simd_lanes(), 0);
#endif
        for (std::size_t active_terms : {0, 1, 7, 8, 9, 16}) {
            polynomial_type outer(2 * block_size, fq12_value_type::zero());
            outer[0] = fq12_value(13);    // A nonzero constant must not move a seven-term sum above the cutoff.
            outer.back() = fq12_value(25);
            // Place active terms at the end of the block, including a lone high-exponent term.
            for (std::size_t j = 1; j <= active_terms; ++j)
                outer[block_size - j] = fq12_value(12 * j + 37);
            BOOST_CHECK_EQUAL(cache.prepared_coefficients().use_for_block(outer, 0), active_terms >= 8);
            BOOST_CHECK(!cache.prepared_coefficients().use_for_block(outer, 1));

            polynomial_type slice(outer.begin(), outer.begin() + block_size), evaluated;
            math::condense(slice);
            math::compose_mod_reference(evaluated, slice, inner, division, context);
            polynomial_type actual(11, fq12_value(61)), expected;
            math::addition(expected, actual, evaluated);
            math::detail::add_composition_block(actual, outer, cache, 0);
            BOOST_CHECK(actual == expected);
            math::subtraction(actual, polynomial_type {}, evaluated);
            math::detail::add_composition_block(actual, outer, cache, 0);
            BOOST_CHECK(actual == polynomial_type {fq12_value_type::zero()});
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_mixes_scalar_blocks_with_dense_batches_and_partial_tails) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    constexpr std::size_t block_size = 17;
    // Five adjacent dense blocks cross a four-block batch boundary. Subsequent blocks alternate between paths,
    // including an all-zero block, a constant block, and seven-term sums immediately below the cutoff.
    const std::array<std::size_t, 14> active_terms = {8, 16, 8, 9, 16, 7, 8, 0, 1, 8, 16, 7, 0, 8};
    polynomial_type divisor(10), inner(9);
    for (std::size_t i = 0; i < divisor.size(); ++i)
        divisor[i] = fq12_value(12 * i + 1);
    for (std::size_t i = 0; i < inner.size(); ++i)
        inner[i] = fq12_value(12 * i + 121);

    for (std::size_t lanes : {0, 4, 8}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = block_size;
        options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
        options.modular_composition_simd_lanes = lanes;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, 8, context);
        const math::polynomial_composition_precomputation<backend_type> cache(inner, block_size * block_size, division,
                                                                              context);
        BOOST_REQUIRE_EQUAL(cache.block_size(), block_size);
        BOOST_REQUIRE(cache.prepared_coefficients().available());
        for (std::size_t tail_size : {1, 2, 8, 9, 17}) {
            polynomial_type outer(active_terms.size() * block_size + tail_size, fq12_value_type::zero());
            for (std::size_t block = 0; block < active_terms.size(); ++block) {
                const std::size_t first = block * block_size;
                if (block != 7)
                    outer[first] = fq12_value(12 * first + 217);
                for (std::size_t j = 1; j <= active_terms[block]; ++j)
                    outer[first + block_size - j] = fq12_value(12 * (first + j) + 229);
                BOOST_CHECK_EQUAL(cache.prepared_coefficients().use_for_block(outer, block), active_terms[block] >= 8);
            }
            for (std::size_t i = active_terms.size() * block_size; i < outer.size(); ++i)
                outer[i] = fq12_value(12 * i + 241);
            BOOST_CHECK_EQUAL(cache.prepared_coefficients().use_for_block(outer, active_terms.size()), tail_size >= 9);
            polynomial_type expected, actual;
            math::compose_mod_reference(expected, outer, inner, division, context);
            math::compose_mod(actual, outer, cache, division, context);
            BOOST_CHECK(actual == expected);
            math::compose_mod(outer, outer, cache, division, context);
            BOOST_CHECK(outer == expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_tiles_match_horner_at_block_and_coefficient_boundaries) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    // Dense full blocks contain eight nonconstant terms, so they still exercise the prepared kernels.
    constexpr std::size_t block_size = 9;
    for (std::size_t degree : {1, 3, 4, 5, 9}) {
        polynomial_type divisor(degree + 1), inner(degree);
        for (std::size_t i = 0; i <= degree; ++i) {
            divisor[i] = fq12_value(i * 12 + 1);
            if (i < degree) {
                inner[i] = fq12_value(i * 12 + 73);
            }
        }
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = block_size;
        options.modular_composition_prepared_byte_limit = 1024 * 1024;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        math::polynomial_divisor_context<backend_type> division(divisor, std::max(std::size_t(1), degree - 1), context);
        const math::polynomial_composition_precomputation<backend_type> cache(inner, 82, division, context);
        BOOST_REQUIRE_EQUAL(cache.block_size(), block_size);
        BOOST_REQUIRE_GT(cache.prepared_coefficient_bytes(), 0);
        for (std::size_t count : {1, 8, 9, 10, 35, 36, 37, 44, 45, 46, 71, 72, 73, 81, 82}) {
            for (bool sparse : {false, true}) {
                polynomial_type outer(count, fq12_value_type::zero());
                for (std::size_t i = 0; i < count; ++i) {
                    if (!sparse || i % block_size == 0 || i + 1 == count) {
                        outer[i] = fq12_value(1 + 12 * i);
                    }
                }
                polynomial_type expected, actual;
                math::compose_mod_reference(expected, outer, inner, division, context);
                math::compose_mod(actual, outer, cache, division, context);
                BOOST_CHECK(actual == expected);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_supports_123_terms_and_the_small_prime_backend) {
    using backend_type = polynomial_arithmetic::bn254_small_prime_fft_backend;
    using polynomial_type = typename backend_type::polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_cached_power_limit = 123;
    options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type(32), options);
    const polynomial_type divisor = {fq12_value(1),  fq12_value(13), fq12_value(25),
                                     fq12_value(37), fq12_value(49), fq12_value(61)};
    const polynomial_type inner = {fq12_value(73), fq12_value(85), fq12_value(97), fq12_value(109), fq12_value(121)};
    math::polynomial_divisor_context<backend_type> division(divisor, 4, context);
    const math::polynomial_composition_precomputation<backend_type> cache(inner, 122 * 122 + 1, division, context);
    BOOST_REQUIRE_EQUAL(cache.block_size(), 123);
    BOOST_REQUIRE_GT(cache.prepared_coefficient_bytes(), 0);
    using reference_backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    polynomial_arithmetic::polynomial_context<reference_backend> reference_context;
    math::polynomial_divisor_context<reference_backend> reference_division(divisor, 4, reference_context);
    for (std::size_t count : {122, 123, 124, 491, 492, 493, 615}) {
        polynomial_type outer(count);
        for (std::size_t i = 0; i < count; ++i) {
            outer[i] = fq12_value(i * 12 + 1);
        }
        polynomial_type expected, actual;
        math::compose_mod_reference(expected, outer, inner, reference_division, reference_context);
        math::compose_mod(actual, outer, cache, division, context);
        BOOST_CHECK(actual == expected);
    }
}

BOOST_AUTO_TEST_CASE(prepared_composition_fallbacks_handle_constant_moduli_and_unsupported_fields) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_prepared_byte_limit = 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    const polynomial_type inner = {fq12_value(1), fq12_value(13)};
    const polynomial_type outer = {fq12_value(25), fq12_value(37), fq12_value(49)};
    math::polynomial_divisor_context<backend_type> constant_division({fq12_value(7)}, 1, context);
    const math::polynomial_composition_precomputation<backend_type> constant_cache(inner, 3, constant_division,
                                                                                   context);
    BOOST_CHECK_EQUAL(constant_cache.prepared_coefficient_bytes(), 0);
    polynomial_type result = outer;
    math::compose_mod(result, result, constant_cache, constant_division, context);
    BOOST_CHECK(result == polynomial_type {fq12_value_type::zero()});

    options.modular_composition_cached_power_limit = 1;
    polynomial_arithmetic::polynomial_context<backend_type> scalar_context(backend_type {}, options);
    math::polynomial_divisor_context<backend_type> division({fq12_value(1), fq12_value(13), fq12_value(25)}, 1,
                                                            scalar_context);
    const math::polynomial_composition_precomputation<backend_type> scalar_cache(inner, 3, division, scalar_context);
    BOOST_CHECK_EQUAL(scalar_cache.prepared_coefficient_bytes(), 0);
    polynomial_type expected;
    math::compose_mod_reference(expected, outer, inner, division, scalar_context);
    math::compose_mod(result, outer, scalar_cache, division, scalar_context);
    BOOST_CHECK(result == expected);

    using base_backend = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    options.modular_composition_cached_power_limit = 3;
    polynomial_arithmetic::polynomial_context<base_backend> base_context(base_backend {}, options);
    math::polynomial_divisor_context<base_backend> base_division({fq_value_type(1), fq_value_type(1)}, 1, base_context);
    const math::polynomial_composition_precomputation<base_backend> base_cache({fq_value_type(2)}, 3, base_division,
                                                                               base_context);
    BOOST_CHECK_EQUAL(base_cache.prepared_coefficient_bytes(), 0);
    base_backend::polynomial_type base_result;
    math::compose_mod(base_result, {fq_value_type(1), fq_value_type(2), fq_value_type(3)}, base_cache, base_division,
                      base_context);
    BOOST_CHECK(base_result == base_backend::polynomial_type {fq_value_type(17)});
}

BOOST_AUTO_TEST_CASE(prepared_composition_handles_degenerate_inputs_and_one_off_aliases) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    polynomial_arithmetic::polynomial_context_options options;
    options.modular_composition_prepared_byte_limit = 1024 * 1024;
    polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
    const polynomial_type divisor = {fq12_value(1), fq12_value(13), fq12_value(25), fq12_value(37), fq12_value(49)};
    math::polynomial_divisor_context<backend_type> division(divisor, 3, context);
    for (const polynomial_type inner : {polynomial_type {fq12_value_type::zero()}, polynomial_type {fq12_value(61)},
                                        polynomial_type {fq12_value(73), fq12_value(85), fq12_value(97)}}) {
        const math::polynomial_composition_precomputation<backend_type> cache(inner, 25, division, context);
        BOOST_REQUIRE_GT(cache.prepared_coefficient_bytes(), 0);
        for (const polynomial_type outer :
             {polynomial_type {fq12_value_type::zero()}, polynomial_type {fq12_value(109)},
              polynomial_type {-inner[0], fq12_value_type::one()},
              polynomial_type {fq12_value(121), fq12_value_type::zero(), fq12_value(133)}}) {
            polynomial_type expected, result;
            math::compose_mod_reference(expected, outer, inner, division, context);
            math::compose_mod(result, outer, cache, division, context);
            BOOST_CHECK(result == expected);
            auto outer_alias = outer;
            math::compose_mod(outer_alias, outer_alias, inner, division, context);
            BOOST_CHECK(outer_alias == expected);
            auto inner_alias = inner;
            math::compose_mod(inner_alias, outer, inner_alias, division, context);
            BOOST_CHECK(inner_alias == expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_simd_panels_match_scalar_at_lane_and_term_boundaries) {
    using polynomial_type = math::polynomial<fq12_value_type>;
    using cache_type = math::detail::prepared_composition_coefficients<fq12_value_type>;
    constexpr std::size_t budget = 16 * 1024 * 1024;
    for (std::size_t width : {3, 4, 5, 7, 8, 9, 15, 17}) {
        std::vector<polynomial_type> powers(18);
        powers[0] = {fq12_value_type::one()};
        for (std::size_t j = 1; j < powers.size(); ++j) {
            powers[j].resize(j % 3 == 0 ? 1 : width);
            for (std::size_t c = 0; c < powers[j].size(); ++c)
                powers[j][c] = fq12_value(12 * (j * width + c) + 1);
        }
        cache_type scalar(powers, budget, 0);
        for (std::size_t lanes : {4, 8}) {
            cache_type candidate(powers, budget, lanes);
            BOOST_REQUIRE(candidate.available());
#if defined(CRYPTO3_HAS_AVX512IFMA)
            BOOST_CHECK_EQUAL(candidate.simd_lanes(), width >= 4 ? lanes : 0);
            if (width >= 4) {
                const std::size_t expected =
                    ((width + lanes - 1) / lanes) * lanes * 23 * 17 * 5 * sizeof(std::uint64_t);
                BOOST_CHECK_EQUAL(candidate.bytes(), expected);
            }
#else
            BOOST_CHECK_EQUAL(candidate.simd_lanes(), 0);
#endif
            BOOST_CHECK_LE(candidate.bytes(), budget);
            const auto bytes = candidate.bytes();
            cache_type exact(powers, bytes, lanes), small(powers, bytes - 1, lanes);
            BOOST_CHECK(exact.available());
            BOOST_CHECK_LE(small.bytes(), bytes - 1);
            auto copied = candidate;
            auto moved = std::move(copied);
            typename cache_type::workspace scratch, reference_scratch;
            std::array<polynomial_type, 4> actual, expected;
            for (bool sparse : {false, true}) {
                polynomial_type outer(5 * powers.size() - 7);
                for (std::size_t i = 0; i < outer.size(); ++i) {
                    outer[i] = sparse && i % 7 ? fq12_value_type::zero() : fq12_value(12 * i + 1);
                }
                // Include a wholly constant block amid active blocks, then a short final block.
                for (std::size_t j = 1; j < powers.size(); ++j)
                    outer[powers.size() + j] = fq12_value_type::zero();
                for (std::size_t first : {0, 1, 4}) {
                    const std::size_t count = std::min(std::size_t(4), std::size_t(5) - first);
                    moved.evaluate_blocks(actual, outer, first, count, scratch);
                    scalar.evaluate_blocks(expected, outer, first, count, reference_scratch);
                    for (std::size_t i = 0; i < count; ++i)
                        BOOST_CHECK(actual[i] == expected[i]);
                }
            }
        }
        BOOST_CHECK_THROW(cache_type(powers, budget, 3), std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(prepared_simd_complete_composition_matches_horner_and_aliases) {
    using backend = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = backend::polynomial_type;
    polynomial_type divisor(10), inner(8), outer(83);
    for (std::size_t i = 0; i < divisor.size(); ++i)
        divisor[i] = fq12_value(12 * i + 1);
    for (std::size_t i = 0; i < inner.size(); ++i)
        inner[i] = fq12_value(12 * i + 121);
    for (std::size_t i = 0; i < outer.size(); ++i)
        outer[i] = i % 3 ? fq12_value(12 * i + 217) : fq12_value_type::zero();
    for (std::size_t lanes : {0, 4, 8}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = 17;
        options.modular_composition_prepared_byte_limit = 4 * 1024 * 1024;
        options.modular_composition_simd_lanes = lanes;
        polynomial_arithmetic::polynomial_context<backend> context(backend {}, options);
        math::polynomial_divisor_context<backend> division(divisor, 8, context);
        const math::polynomial_composition_precomputation<backend> cache(inner, 257, division, context);
        polynomial_type expected, actual = outer;
        math::compose_mod_reference(expected, outer, inner, division, context);
        math::compose_mod(actual, actual, cache, division, context);
        BOOST_CHECK(actual == expected);
    }
}

BOOST_AUTO_TEST_CASE(prepared_simd_configuration_keeps_six_limb_fields_on_scalar_kernel) {
    using value = fields::fp12_2over3over2<fields::bls12<381>>::value_type;
    using polynomial = math::polynomial<value>;
    using cache_type = math::detail::prepared_composition_coefficients<value>;
    std::vector<polynomial> powers(9, polynomial(4, value::one()));
    powers[0] = {value::one()};
    const cache_type cache(powers, 1 << 20, 8);
    BOOST_REQUIRE(cache.available());
    BOOST_CHECK_EQUAL(cache.simd_lanes(), 0);
    typename cache_type::workspace scratch;
    std::array<polynomial, 4> output;
    cache.evaluate_blocks(output, polynomial(9, value::one()), 0, 1, scratch);
    BOOST_CHECK(output[0] == polynomial({value(9), value(8), value(8), value(8)}));
}

BOOST_AUTO_TEST_SUITE_END()
