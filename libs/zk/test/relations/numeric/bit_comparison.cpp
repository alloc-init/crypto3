//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init
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

#define BOOST_TEST_MODULE bit_comparison_test

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/zk/snark/arithmetization/bit_comparison.hpp>

namespace {
    using boost::multiprecision::cpp_int;
    using nil::crypto3::zk::snark::for_each_bit_comparison;

    // Evaluate the emitted equations modulo a small prime, independently of circuit/witness containers.
    // next_marker can supply adversarial values instead of the honest product.
    template<typename Integer, typename NextMarker>
    bool comparison_satisfied(const Integer &value, std::size_t bit_count, const Integer &bound,
                              unsigned characteristic, NextMarker &&next_marker, bool square_markers = false) {
        unsigned marker = 1;
        bool satisfied = true;
        for_each_bit_comparison(
            bound, bit_count, characteristic,
            [&](std::size_t bit) { marker = boost::multiprecision::bit_test(value, bit); },
            [&](std::size_t bit) {
                const unsigned product = marker * boost::multiprecision::bit_test(value, bit) % characteristic;
                marker = next_marker(product);
                const unsigned output = square_markers ? marker * marker % characteristic : marker;
                satisfied &= output == product;
            },
            [&](std::size_t begin, std::size_t end) {
                unsigned sum = 0;
                for (std::size_t bit = begin; bit < end; ++bit) {
                    sum = (sum + boost::multiprecision::bit_test(value, bit)) % characteristic;
                }
                satisfied &= marker * sum % characteristic == 0;
            });
        return satisfied;
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(bit_comparison_test_suite)

BOOST_AUTO_TEST_CASE(exhaustive_small_integers_and_widths) {
    // Includes empty inputs, vacuous bounds, leading/trailing zeros and powers of two.
    // A small characteristic also exercises zero-run splitting; the bound can exceed it.
    for (std::size_t width = 0; width <= 8; ++width) {
        for (unsigned bound = 1; bound <= 256; ++bound) {
            for (unsigned value = 0; value < (1u << width); ++value) {
                BOOST_TEST_CONTEXT("width=" << width << ", bound=" << bound << ", value=" << value) {
                    BOOST_CHECK_EQUAL(
                        comparison_satisfied(value, width, bound, 5, [](unsigned product) { return product; }),
                        value < bound);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(bn254_counts_are_constexpr_and_bound_is_independent_of_field) {
    using scalar_field = nil::crypto3::algebra::fields::alt_bn128_scalar_field<254>;
    using base_field = nil::crypto3::algebra::fields::alt_bn128_base_field<254>;
    constexpr auto scalar_counts = for_each_bit_comparison(
        scalar_field::modulus, 254, scalar_field::modulus, [](std::size_t) { }, [](std::size_t) { },
        [](std::size_t, std::size_t) { });
    constexpr auto scalar_in_base_counts = for_each_bit_comparison(
        scalar_field::modulus, 254, base_field::modulus, [](std::size_t) { }, [](std::size_t) { },
        [](std::size_t, std::size_t) { });
    constexpr auto base_counts = for_each_bit_comparison(
        base_field::modulus, 254, base_field::modulus, [](std::size_t) { }, [](std::size_t) { },
        [](std::size_t, std::size_t) { });
    BOOST_CHECK_EQUAL(scalar_counts.marker_count, 99);
    BOOST_CHECK_EQUAL(scalar_counts.product_count, 153);
    BOOST_CHECK_EQUAL(scalar_in_base_counts.marker_count, 99);
    BOOST_CHECK_EQUAL(scalar_in_base_counts.product_count, 153);
    BOOST_CHECK_EQUAL(base_counts.marker_count, 109);
    BOOST_CHECK_EQUAL(base_counts.product_count, 171);
}

BOOST_AUTO_TEST_CASE(callback_order_and_bit_positions) {
    std::vector<std::string> operations;
    // 52 = 00110100 in eight bits. The final comparison is at bit 2, its lowest one.
    const auto counts = for_each_bit_comparison(
        52u, 8, 3u, [&](std::size_t bit) { operations.push_back("reuse " + std::to_string(bit)); },
        [&](std::size_t bit) { operations.push_back("update " + std::to_string(bit)); },
        [&](std::size_t begin, std::size_t end) {
            operations.push_back("zero " + std::to_string(begin) + " " + std::to_string(end));
        });
    const std::vector<std::string> expected = {"zero 6 8", "reuse 5", "update 4", "zero 3 4", "zero 2 3"};
    BOOST_CHECK(operations == expected);
    BOOST_CHECK_EQUAL(counts.marker_count, 1);
    BOOST_CHECK_EQUAL(counts.product_count, 4);
}

BOOST_AUTO_TEST_CASE(alternative_marker_values_preserve_the_strict_bound) {
    // Bound 7 has one fresh marker. Test every possible field value with either marker equation.
    for (const bool square_markers : {false, true}) {
        for (unsigned value = 0; value < 8; ++value) {
            bool accepted = false;
            for (unsigned marker = 0; marker < 7; ++marker) {
                BOOST_TEST_CONTEXT("squared=" << square_markers << ", value=" << value << ", marker=" << marker) {
                    const bool valid =
                        comparison_satisfied(value, 3, 7u, 7, [=](unsigned) { return marker; }, square_markers);
                    BOOST_CHECK(!valid || value < 7);
                    accepted |= valid;
                }
            }
            BOOST_CHECK_EQUAL(accepted, value < 7);
        }
    }

    // Over F17, (-1)^2 = 1 and 4^2 = -1. Squared markers need not be Boolean.
    const std::vector<unsigned> markers = {16, 4};
    for (const unsigned value : {14, 15}) {
        std::size_t next = 0;
        BOOST_CHECK_EQUAL(
            comparison_satisfied(value, 4, 15u, 17, [&](unsigned) { return markers.at(next++); }, true), value < 15);
    }
    std::size_t next = 0;
    BOOST_CHECK(!comparison_satisfied(14u, 4, 15u, 17, [&](unsigned) { return markers.at(next++); }));
}

BOOST_AUTO_TEST_CASE(zero_runs_do_not_wrap_modulo_the_characteristic) {
    // The 251 zeros in the bound must be split: 251 Boolean ones would sum to zero over F251.
    const cpp_int bound = (cpp_int(1) << 252) + 1;
    const auto counts = for_each_bit_comparison(
        bound, 253, 251u, [](std::size_t) { }, [](std::size_t) { }, [](std::size_t, std::size_t) { });
    BOOST_CHECK_EQUAL(counts.marker_count, 0);
    BOOST_CHECK_EQUAL(counts.product_count, 3);
    const std::vector<cpp_int> values = {0, bound - 1, bound, (cpp_int(1) << 253) - 2};
    for (const auto &value : values) {
        BOOST_CHECK_EQUAL(comparison_satisfied(value, 253, bound, 251, [](unsigned product) { return product; }),
                          value < bound);
    }

    // In characteristic two, even a run of two zeros must be split into single-bit checks.
    for (const unsigned value : {0, 8, 9, 14, 15}) {
        BOOST_CHECK_EQUAL(comparison_satisfied(value, 4, 9u, 2, [](unsigned product) { return product; }), value < 9);
    }
}

BOOST_AUTO_TEST_CASE(invalid_parameters_fail_before_callbacks_and_vacuous_bounds_emit_nothing) {
    std::size_t calls = 0;
    const auto record_call = [&](auto...) { ++calls; };
    for (const std::size_t width : {0, 8}) {
        for (const int bound : {0, -1}) {
            BOOST_CHECK_THROW(for_each_bit_comparison(bound, width, 5, record_call, record_call, record_call),
                              std::invalid_argument);
        }
        for (const int characteristic : {-1, 0, 1}) {
            BOOST_CHECK_THROW(
                for_each_bit_comparison(256u, width, characteristic, record_call, record_call, record_call),
                std::invalid_argument);
        }
    }
    BOOST_CHECK_EQUAL(calls, 0);
    for (const std::size_t width : {0, 3, 8}) {
        const auto counts = for_each_bit_comparison(256u, width, 2u, record_call, record_call, record_call);
        BOOST_CHECK_EQUAL(counts.marker_count, 0);
        BOOST_CHECK_EQUAL(counts.product_count, 0);
    }
    BOOST_CHECK_EQUAL(calls, 0);
}

BOOST_AUTO_TEST_SUITE_END()
