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

#include <array>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/random/mersenne_twister.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/math/algorithms/mixed_radix_fft.hpp>
#include <nil/crypto3/math/polynomial/operations/basic_operations.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>

namespace {

    using clock_type = std::chrono::steady_clock;
    using bn254_fq = nil::crypto3::algebra::fields::alt_bn128_base_field<254>;
    using bn254_fq12 = nil::crypto3::algebra::fields::fp12_2over3over2<nil::crypto3::algebra::fields::alt_bn128<254>>;
    using fq12_value_type = bn254_fq12::value_type;
    using backend_type = nil::crypto3::math::polynomial_arithmetic::mixed_radix_backend<bn254_fq, fq12_value_type>;
    using context_type = nil::crypto3::math::polynomial_arithmetic::polynomial_context<backend_type>;
    using polynomial_type = typename context_type::polynomial_type;

    constexpr std::size_t odd_smooth_order = 3 * 3 * 29 * 67;
    constexpr std::size_t even_smooth_order = 2 * odd_smooth_order;

    template<typename Function>
    double elapsed_milliseconds(Function &&function) {
        const auto start = clock_type::now();
        function();
        const auto finish = clock_type::now();
        return std::chrono::duration<double, std::milli>(finish - start).count();
    }

    std::vector<fq12_value_type> random_values(std::size_t size, boost::random::mt19937 &rng) {
        std::vector<fq12_value_type> values;
        values.reserve(size);
        for (std::size_t i = 0; i < size; ++i) {
            values.push_back(nil::crypto3::algebra::random_element<bn254_fq12>(rng));
        }
        return values;
    }

    void run_benchmark(std::size_t size, boost::random::mt19937 &rng) {
        const auto plan_start = clock_type::now();
        const nil::crypto3::math::mixed_radix_fft_plan<bn254_fq> plan(size);
        const auto plan_finish = clock_type::now();
        const double plan_ms = std::chrono::duration<double, std::milli>(plan_finish - plan_start).count();

        const auto backend_start = clock_type::now();
        context_type context {backend_type(size)};
        const auto backend_finish = clock_type::now();
        const double backend_ms = std::chrono::duration<double, std::milli>(backend_finish - backend_start).count();

        const std::vector<fq12_value_type> coefficients = random_values(size, rng);
        std::vector<fq12_value_type> transformed = coefficients;
        const double forward_ms = elapsed_milliseconds([&]() { plan.fft(transformed); });
        const double inverse_ms = elapsed_milliseconds([&]() { plan.inverse_fft(transformed); });
        if (transformed != coefficients) {
            throw std::runtime_error("mixed-radix FFT round trip failed");
        }

        const std::size_t left_size = (size + 1) / 2;
        const std::size_t right_size = size - left_size + 1;
        const polynomial_type left(random_values(left_size, rng));
        const polynomial_type right(random_values(right_size, rng));
        polynomial_type product;
        const double multiplication_ms =
            elapsed_milliseconds([&]() { nil::crypto3::math::multiplication(product, left, right, context); });
        if (product.size() != size) {
            throw std::runtime_error("mixed-radix multiplication returned the wrong size");
        }

        std::cout << size << ',' << plan_ms << ',' << backend_ms << ',' << forward_ms << ',' << inverse_ms << ','
                  << multiplication_ms << '\n';
    }

}    // namespace

int main() {
    boost::random::mt19937 rng(0xC0FFEE);

    std::cout << "Mixed-radix Fq12 benchmark; milliseconds; no pass/fail thresholds\n";
    std::cout << "size,fft_plan,backend,forward,inverse,multiplication\n";
    std::cout << std::fixed << std::setprecision(3);

    for (const std::size_t size : std::array<std::size_t, 2> {odd_smooth_order, even_smooth_order}) {
        run_benchmark(size, rng);
    }
}
