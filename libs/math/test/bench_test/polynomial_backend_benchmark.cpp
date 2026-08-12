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

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <boost/random/mersenne_twister.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/math/polynomial/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/schoolbook_backend.hpp>

namespace {

    namespace fields = nil::crypto3::algebra::fields;
    namespace polynomial_arithmetic = nil::crypto3::math::polynomial_arithmetic;

    using clock_type = std::chrono::steady_clock;
    using fq_field_type = fields::alt_bn128_base_field<254>;
    using fq12_field_type = fields::fp12_2over3over2<fields::alt_bn128<254>>;
    using value_type = fq12_field_type::value_type;
    using coefficient_vector = polynomial_arithmetic::coefficient_vector<value_type>;

    constexpr std::chrono::milliseconds minimum_measurement_time(100);
    constexpr std::size_t maximum_repetitions = 1024;

    // These nontrivial divisors of 2 * 3^2 * 29 * 67 cover the radix combinations
    // near the expected schoolbook/FFT crossover without making quadratic
    // measurements at the largest supported orders.
    constexpr std::array<std::size_t, 15> transform_sizes = {
        3, 6, 9, 18, 29, 58, 67, 87, 134, 174, 201, 261, 402, 522, 603,
    };

    coefficient_vector random_values(std::size_t size, boost::random::mt19937 &rng) {
        coefficient_vector values;
        values.reserve(size);
        for (std::size_t i = 0; i < size; ++i) {
            values.push_back(nil::crypto3::algebra::random_element<fq12_field_type>(rng));
        }
        return values;
    }

    template<typename Function>
    double average_milliseconds(Function &&function) {
        std::size_t repetitions = 1;
        while (true) {
            const auto start = clock_type::now();
            for (std::size_t i = 0; i < repetitions; ++i) {
                function();
            }
            const auto finish = clock_type::now();
            const auto elapsed = finish - start;
            if (elapsed >= minimum_measurement_time || repetitions == maximum_repetitions) {
                return std::chrono::duration<double, std::milli>(elapsed).count() / repetitions;
            }
            repetitions = std::min(2 * repetitions, maximum_repetitions);
        }
    }

    template<typename SchoolbookOperation, typename MixedRadixOperation>
    void benchmark_operation(std::size_t transform_size, std::string_view operation,
                             SchoolbookOperation &&schoolbook_operation, MixedRadixOperation &&mixed_radix_operation,
                             coefficient_vector &schoolbook_output, coefficient_vector &mixed_radix_output) {
        schoolbook_operation();
        mixed_radix_operation();
        if (schoolbook_output != mixed_radix_output) {
            throw std::runtime_error("polynomial backends returned different results");
        }

        const double schoolbook_ms = average_milliseconds(schoolbook_operation);
        const double mixed_radix_ms = average_milliseconds(mixed_radix_operation);
        std::cout << transform_size << ',' << operation << ',' << schoolbook_ms << ',' << mixed_radix_ms << ','
                  << schoolbook_ms / mixed_radix_ms << '\n';
    }

    void run_benchmark(std::size_t transform_size, boost::random::mt19937 &rng) {
        const std::size_t left_size = (transform_size + 1) / 2;
        const std::size_t right_size = transform_size - left_size + 1;
        const std::size_t low_coefficient_count = (transform_size + 1) / 2;
        const coefficient_vector left = random_values(left_size, rng);
        const coefficient_vector right = random_values(right_size, rng);

        polynomial_arithmetic::schoolbook_backend<value_type> schoolbook;
        polynomial_arithmetic::mixed_radix_backend<fq_field_type, value_type> mixed_radix(transform_size);
        coefficient_vector schoolbook_output;
        coefficient_vector mixed_radix_output;

        benchmark_operation(
            transform_size, "multiply", [&]() { schoolbook.multiply(schoolbook_output, left, right); },
            [&]() { mixed_radix.multiply(mixed_radix_output, left, right); }, schoolbook_output, mixed_radix_output);
        benchmark_operation(
            transform_size, "square", [&]() { schoolbook.square(schoolbook_output, left); },
            [&]() { mixed_radix.square(mixed_radix_output, left); }, schoolbook_output, mixed_radix_output);
        benchmark_operation(
            transform_size, "multiply_low",
            [&]() { schoolbook.multiply_low(schoolbook_output, left, right, low_coefficient_count); },
            [&]() { mixed_radix.multiply_low(mixed_radix_output, left, right, low_coefficient_count); },
            schoolbook_output, mixed_radix_output);
    }

}    // namespace

int main() {
    boost::random::mt19937 rng(0xC0FFEE);

    std::cout << "BN254 Fq12 polynomial-backend benchmark; milliseconds per operation; no pass/fail thresholds\n";
    std::cout << "transform_size,operation,schoolbook,mixed_radix,schoolbook_over_mixed_radix\n";
    std::cout << std::fixed << std::setprecision(3);

    for (const std::size_t transform_size : transform_sizes) {
        run_benchmark(transform_size, rng);
    }
}
