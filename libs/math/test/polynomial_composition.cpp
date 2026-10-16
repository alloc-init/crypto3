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

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_composition.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>

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

BOOST_AUTO_TEST_CASE(coefficient_blocks_preserve_ragged_powers_and_accumulator_tails) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    polynomial_arithmetic::polynomial_context<backend_type> context;
    polynomial_type divisor(10, fq_value_type::zero());
    divisor.front() = divisor.back() = fq_value_type::one();
    const math::polynomial_divisor_context<backend_type> modulus(divisor, 8, context);
    const polynomial_type inner {fq_value_type::zero(), fq_value_type::one()};
    const math::polynomial_composition_precomputation<backend_type> cache(inner, 18, modulus, context);
    BOOST_REQUIRE_EQUAL(cache.block_size(), 5);
    for (std::size_t j = 0; j < cache.block_size(); ++j) {
        BOOST_CHECK_EQUAL(cache.baby_step(j).size(), j + 1);
    }

    const polynomial_type initial(8, fq_value_type(17));
    for (std::size_t size : {1, 4, 5, 6, 18}) {
        polynomial_type outer(size, fq_value_type::zero());
        for (std::size_t i = 0; i < size; ++i) {
            outer[i] = i % 3 == 0 && i + 1 < size ? fq_value_type::zero() : fq_value_type(i + 1);
        }
        const auto blocks = 1 + (size - 1) / cache.block_size();
        for (std::size_t block = blocks; block-- > 0;) {
            // With inner = X and no wraparound, each block is exactly its coefficient slice.
            auto expected = initial;
            const auto first = block * cache.block_size();
            for (std::size_t j = 0; j < std::min(cache.block_size(), size - first); ++j) {
                expected[j] += outer[first + j];
            }
            auto actual = initial;
            math::detail::add_composition_block(actual, outer, cache, block);
            BOOST_CHECK(actual == expected);
            BOOST_CHECK_EQUAL(actual.size(), initial.size());
        }
    }
}

BOOST_AUTO_TEST_CASE(coefficient_blocks_match_fp12_horner_for_sparse_tails_and_cancellation) {
    using backend_type = polynomial_arithmetic::schoolbook_backend<fq12_value_type>;
    using polynomial_type = typename backend_type::polynomial_type;
    const polynomial_type divisor {fq12_value(1),  fq12_value(13), fq12_value(25), fq12_value(37),
                                   fq12_value(49), fq12_value(61), fq12_value(73), fq12_value(85)};
    const polynomial_type dense_outer {fq12_value(97),  fq12_value(109), fq12_value(121), fq12_value(133),
                                       fq12_value(145), fq12_value(157), fq12_value(169), fq12_value(181),
                                       fq12_value(193), fq12_value(205), fq12_value(217), fq12_value(229),
                                       fq12_value(241)};
    const polynomial_type sparse_outer {fq12_value(253), fq12_value_type::zero(), fq12_value_type::zero(),
                                        fq12_value_type::zero(), fq12_value(265)};
    const polynomial_type initial {fq12_value(277), fq12_value(289), fq12_value(301), fq12_value(313)};
    for (std::size_t limit : {1, 2, 4}) {
        polynomial_arithmetic::polynomial_context_options options;
        options.modular_composition_cached_power_limit = limit;
        polynomial_arithmetic::polynomial_context<backend_type> context(backend_type {}, options);
        const math::polynomial_divisor_context<backend_type> modulus(divisor, 6, context);
        for (const auto &inner : {polynomial_type {fq12_value(325)}, polynomial_type {fq12_value(337), fq12_value(349)},
                                  polynomial_type {fq12_value(361), fq12_value(373), fq12_value(385), fq12_value(397),
                                                   fq12_value(409), fq12_value(421), fq12_value(433)}}) {
            const math::polynomial_composition_precomputation<backend_type> cache(inner, dense_outer.size(), modulus,
                                                                                  context);
            BOOST_REQUIRE_EQUAL(cache.block_size(), limit);
            for (const auto &outer : {polynomial_type {fq12_value_type::zero()}, polynomial_type {fq12_value(445)},
                                      sparse_outer, dense_outer}) {
                const auto blocks = 1 + (outer.size() - 1) / cache.block_size();
                for (std::size_t block = blocks; block-- > 0;) {
                    const auto first = block * cache.block_size();
                    const auto end = std::min(first + cache.block_size(), outer.size());
                    polynomial_type block_polynomial(outer.begin() + first, outer.begin() + end);
                    math::condense(block_polynomial);
                    polynomial_type block_value;
                    math::compose_mod_reference(block_value, block_polynomial, inner, modulus, context);
                    for (const auto &accumulator : {polynomial_type {}, initial}) {
                        polynomial_type expected;
                        math::addition(expected, accumulator, block_value);
                        auto actual = accumulator;
                        math::detail::add_composition_block(actual, outer, cache, block);
                        BOOST_CHECK(actual == expected);
                    }
                    auto cancelled = block_value;
                    for (auto &coefficient : cancelled) {
                        coefficient = -coefficient;
                    }
                    math::detail::add_composition_block(cancelled, outer, cache, block);
                    BOOST_CHECK_EQUAL(cancelled.size(), 1);
                    BOOST_CHECK(cancelled.front().is_zero());
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
