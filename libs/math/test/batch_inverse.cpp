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

#define BOOST_TEST_MODULE batch_inverse_test

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/math/algorithms/batch_inverse.hpp>

using namespace nil::crypto3;

namespace {

    using fq_field = algebra::fields::alt_bn128<254>;
    using fq12_field = algebra::fields::fp12_2over3over2<fq_field>;

    template<typename FieldType>
    void check_empty_input() {
        using value_type = typename FieldType::value_type;

        const std::vector<value_type> values;
        BOOST_CHECK(math::batch_inverse_nonzero(values).empty());
    }

    template<typename FieldType>
    void check_singleton_input() {
        using value_type = typename FieldType::value_type;

        const value_type value = value_type::one().doubled() + value_type::one();
        const std::vector<value_type> values = {value};
        const auto inverses = math::batch_inverse_nonzero(values);

        BOOST_REQUIRE_EQUAL(inverses.size(), 1);
        BOOST_CHECK_EQUAL(inverses[0], values[0].inversed());
    }

    template<typename FieldType>
    void check_random_input() {
        using value_type = typename FieldType::value_type;

        boost::random::mt19937 rng(0xB47C4u);
        std::vector<value_type> values(257);
        for (value_type &value : values) {
            do {
                value = algebra::random_element<FieldType>(rng);
            } while (value.is_zero());
        }
        const auto original_values = values;

        const auto inverses = math::batch_inverse_nonzero(values);

        BOOST_REQUIRE_EQUAL(inverses.size(), values.size());
        BOOST_CHECK(values == original_values);
        for (std::size_t i = 0; i < values.size(); ++i) {
            BOOST_CHECK_EQUAL(values[i] * inverses[i], value_type::one());
        }
    }

    template<typename FieldType>
    void check_zero_inputs_are_rejected() {
        using value_type = typename FieldType::value_type;

        for (const std::size_t zero_index : {std::size_t(0), std::size_t(2), std::size_t(4)}) {
            const value_type nonzero = value_type::one().doubled() + value_type::one();
            std::vector<value_type> values(5, nonzero);
            values[zero_index] = value_type::zero();
            BOOST_CHECK_THROW(static_cast<void>(math::batch_inverse_nonzero(values)), std::invalid_argument);
        }
    }

    class counted_value {
    public:
        using value_type = fq_field::value_type;

        static std::size_t multiplication_count;
        static std::size_t inversion_count;

        counted_value() = default;
        counted_value(const value_type &value) : value(value) {
        }
        counted_value(unsigned value) : value(value) {
        }

        bool is_zero() const {
            return value.is_zero();
        }

        counted_value inversed() const {
            ++inversion_count;
            return value.inversed();
        }

        counted_value operator*(const counted_value &other) const {
            ++multiplication_count;
            return value * other.value;
        }

        static void reset_counts() {
            multiplication_count = 0;
            inversion_count = 0;
        }

    private:
        value_type value;
    };

    std::size_t counted_value::multiplication_count = 0;
    std::size_t counted_value::inversion_count = 0;

}    // namespace

BOOST_AUTO_TEST_SUITE(batch_inverse_test_suite)

BOOST_AUTO_TEST_CASE(empty_input) {
    check_empty_input<fq_field>();
    check_empty_input<fq12_field>();
}

BOOST_AUTO_TEST_CASE(singleton_input) {
    check_singleton_input<fq_field>();
    check_singleton_input<fq12_field>();
}

BOOST_AUTO_TEST_CASE(random_input) {
    check_random_input<fq_field>();
    check_random_input<fq12_field>();
}

BOOST_AUTO_TEST_CASE(zero_inputs_are_rejected) {
    check_zero_inputs_are_rejected<fq_field>();
    check_zero_inputs_are_rejected<fq12_field>();
}

BOOST_AUTO_TEST_CASE(uses_one_inversion_and_optimal_multiplication_count) {
    constexpr std::size_t input_size = 17;
    std::vector<counted_value> values;
    values.reserve(input_size);
    for (std::size_t i = 0; i < input_size; ++i) {
        values.emplace_back(static_cast<unsigned>(i + 1));
    }

    counted_value::reset_counts();
    const auto inverses = math::batch_inverse_nonzero(values);

    BOOST_REQUIRE_EQUAL(inverses.size(), input_size);
    BOOST_CHECK_EQUAL(counted_value::inversion_count, 1);
    BOOST_CHECK_EQUAL(counted_value::multiplication_count, 3 * (input_size - 1));
}

BOOST_AUTO_TEST_SUITE_END()
