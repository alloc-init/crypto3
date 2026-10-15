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
#include <bit>
#include <charconv>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include <boost/random/mersenne_twister.hpp>
#include <boost/version.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/alt_bn128.hpp>
#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/mixed_radix_backend.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/operations/evaluate.hpp>
#include <nil/crypto3/math/polynomial/quotient_ring/polynomial_modular_arithmetic.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace pa = math::polynomial_arithmetic;
    namespace fs = std::filesystem;
    using small_backend = pa::bn254_small_prime_fft_backend;
    using field_type = small_backend::field_type;
    using value_type = small_backend::value_type;
    using polynomial = small_backend::polynomial_type;
    using root_field = nil::crypto3::algebra::fields::alt_bn128_base_field<254>;
    using native_backend = pa::mixed_radix_backend<root_field, value_type>;
    constexpr std::size_t native_max_order = 34974;
    constexpr unsigned seed = 0xC0FFEE;

    void require(bool condition, const char *message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    struct timing {
        double wall_ms, cpu_ms;
    };

    template<typename Function>
    timing measure(Function &&function) {
        const auto cpu_start = std::clock();
        const auto start = std::chrono::steady_clock::now();
        function();
        const auto finish = std::chrono::steady_clock::now();
        const auto cpu_finish = std::clock();
        require(cpu_start != std::clock_t(-1) && cpu_finish != std::clock_t(-1), "CPU clock unavailable");
        return {std::chrono::duration<double, std::milli>(finish - start).count(),
                1000.0 * (cpu_finish - cpu_start) / CLOCKS_PER_SEC};
    }

    std::size_t native_order(std::size_t required) {
        for (std::size_t order = required; order <= native_max_order; ++order) {
            if (native_max_order % order == 0) {
                return order;
            }
        }
        throw std::invalid_argument("no supported native transform order");
    }

    struct options {
        std::size_t degree = 14938, trials = 5, warmups = 2;
        bool check = false;
        fs::path input_directory;
    };

    // Export canonical coordinates outside timed regions so the runner can hash the exact inputs.
    void write_input(const polynomial &input, const fs::path &path) {
        std::ofstream output;
        output.exceptions(std::ios::failbit | std::ios::badbit);
        output.open(path, std::ios::out | std::ios::noreplace);
        output << input.size() << '\n' << std::hex;
        for (const auto &coefficient : input) {
            for (std::size_t c = 0; c < field_type::arity; ++c) {
                output << coefficient.data[c / 6].data[(c % 6) / 2].data[c % 2].data << '\n';
            }
        }
    }

    struct fixture {
        polynomial h, a, b, product, square, low, quotient, remainder, square_remainder;
        std::string name;

        fixture(const options &settings, bool sparse) :
            name(sparse ? "synthetic_nonmonic_sparse" : "synthetic_monic_dense") {
            boost::random::mt19937 rng(seed);
            const auto generate = [&](std::size_t size) {
                polynomial result(size, value_type::zero());
                for (auto &coefficient : result) {
                    coefficient = nil::crypto3::algebra::random_element<field_type>(rng);
                }
                require(!result.back().is_zero(), "zero leading coefficient in deterministic fixture");
                return result;
            };
            h = generate(settings.degree + 1);
            a = generate(settings.degree);
            b = generate(settings.degree);
            if (sparse) {
                for (std::size_t i = 1; i + 1 < a.size(); ++i) {
                    a[i] = value_type::zero();
                }
            } else {
                h.back() = value_type::one();
            }
            if (!settings.input_directory.empty()) {
                fs::create_directories(settings.input_directory);
                const auto prefix = std::to_string(settings.degree) + "-" + name + "-";
                write_input(h, settings.input_directory / (prefix + "h.txt"));
                write_input(a, settings.input_directory / (prefix + "a.txt"));
                write_input(b, settings.input_directory / (prefix + "b.txt"));
            }
        }

        template<typename Backend>
        void reference(pa::polynomial_context<Backend> &context, bool classical) {
            context.multiply(product, a, b);
            context.square(square, a);
            polynomial square_quotient;
            if (classical) {
                math::division(quotient, remainder, product, h);
                math::division(square_quotient, square_remainder, square, h);
            } else {
                math::polynomial_divisor_context<Backend> divisor(h, h.size() - 2, context);
                math::divrem(quotient, remainder, product, divisor, context);
                math::divrem(square_quotient, square_remainder, square, divisor, context);
            }
            low = product;
            low.resize(a.size());
            math::condense(low);
            polynomial reconstructed;
            context.multiply(reconstructed, quotient, h);
            math::addition(reconstructed, reconstructed, remainder);
            require(reconstructed == product && remainder.size() < h.size(), "product division identity failed");
            context.multiply(reconstructed, square_quotient, h);
            math::addition(reconstructed, reconstructed, square_remainder);
            require(reconstructed == square && square_remainder.size() < h.size(), "square division identity failed");
            for (unsigned point : {0u, 1u, 7u}) {
                const value_type x(point);
                const auto av = math::evaluate_polynomial(a, x, a.size());
                const auto bv = math::evaluate_polynomial(b, x, b.size());
                require(math::evaluate_polynomial(product, x, product.size()) == av * bv,
                        "product evaluation identity failed");
                require(math::evaluate_polynomial(square, x, square.size()) == av * av,
                        "square evaluation identity failed");
            }
        }
    };

    enum class operation { multiply, square, low, prepared_low, cyclic, divrem, mulmod, squaremod };
    constexpr std::array<std::string_view, 8> operation_names {
        "multiply", "square", "multiply_low", "prepared_low", "cyclic_remainder", "divrem", "mulmod", "squaremod"};

    template<typename Backend>
    class benchmark_case {
        using context_type = pa::polynomial_context<Backend>;
        using divisor_type = math::polynomial_divisor_context<Backend>;
        const fixture &data_;
        std::string_view name_;
        std::size_t order_;
        std::optional<context_type> context_;
        std::optional<divisor_type> divisor_;
        std::optional<typename Backend::prepared_low_product_type> low_;
        std::optional<typename Backend::prepared_cyclic_remainder_type> cyclic_;
        polynomial output_, quotient_;

        void print(std::string_view phase, std::size_t trial, std::string_view operation, timing elapsed,
                   std::size_t transform = 0, std::size_t bytes = 0) const {
            std::cout << name_ << ',' << data_.a.size() << ',' << data_.name << ',' << order_ << ',' << phase << ','
                      << trial << ',' << operation << ',' << elapsed.wall_ms << ',' << elapsed.cpu_ms << ','
                      << transform << ',' << bytes << '\n';
        }

    public:
        benchmark_case(const fixture &data, std::string_view name, std::size_t order) :
            data_(data), name_(name), order_(order) {
            print("setup", 0, "backend", measure([&] { context_.emplace(Backend(order)); }));
            print("setup", 0, "divisor_context",
                  measure([&] { divisor_.emplace(data.h, std::max(std::size_t(1), data.a.size() - 1), *context_); }));
            const auto low_time =
                measure([&] { low_ = context_->prepare_low_product(data.b, data.a.size(), data.a.size()); });
            require(low_.has_value(), "fixed low-product preparation failed");
            print("setup", 0, "low_operand", low_time, low_->transform_size(), low_->storage_bytes());
            const auto cyclic_time =
                measure([&] { cyclic_ = context_->prepare_cyclic_remainder(data.h, data.a.size() - 1); });
            print("setup", 0, "cyclic_divisor", cyclic_time, cyclic_ ? cyclic_->transform_size() : 0,
                  cyclic_ ? cyclic_->storage_bytes() : 0);
        }

        void sample(operation op, std::string_view phase, std::size_t trial) {
            const auto label = operation_names[static_cast<std::size_t>(op)];
            if (op == operation::cyclic && !cyclic_) {
                if (trial == 0) {
                    print("unsupported", 0, label, {0, 0});
                }
                return;
            }
            bool accepted = true;
            const auto elapsed = measure([&] {
                switch (op) {
                    case operation::multiply:
                        context_->multiply(output_, data_.a, data_.b);
                        break;
                    case operation::square:
                        context_->square(output_, data_.a);
                        break;
                    case operation::low:
                        context_->multiply_low(output_, data_.a, data_.b, data_.a.size());
                        break;
                    case operation::prepared_low:
                        accepted = context_->try_multiply_low_prepared(output_, data_.a, *low_, data_.a.size());
                        break;
                    case operation::cyclic:
                        accepted = context_->try_cyclic_remainder(output_, data_.product, data_.quotient, *cyclic_);
                        break;
                    case operation::divrem:
                        math::divrem(quotient_, output_, data_.product, *divisor_, *context_);
                        break;
                    case operation::mulmod:
                        math::mulmod(output_, data_.a, data_.b, *divisor_, *context_);
                        break;
                    case operation::squaremod:
                        math::squaremod(output_, data_.a, *divisor_, *context_);
                        break;
                }
            });
            // Consume and check complete outputs outside the timer on every sample, including warmups.
            require(accepted, "prepared operation declined");
            const polynomial &expected = op == operation::multiply                             ? data_.product :
                                         op == operation::square                               ? data_.square :
                                         op == operation::low || op == operation::prepared_low ? data_.low :
                                         op == operation::squaremod ? data_.square_remainder :
                                                                      data_.remainder;
            require(output_ == expected && !output_.empty() && (output_.size() == 1 || !output_.back().is_zero()),
                    "benchmark result mismatch");
            if (op == operation::divrem) {
                require(quotient_ == data_.quotient, "benchmark quotient mismatch");
            }
            print(phase, trial, label, elapsed,
                  op == operation::prepared_low ? low_->transform_size() :
                  op == operation::cyclic       ? cyclic_->transform_size() :
                                                  0);
        }
    };

    void run(const options &settings, bool sparse) {
        fixture data(settings, sparse);
        const auto order = native_order(2 * settings.degree - 1);
        if (settings.degree <= 67) {
            pa::polynomial_context<pa::schoolbook_backend<value_type>> reference;
            data.reference(reference, true);
        } else {
            pa::polynomial_context<native_backend> reference {native_backend(order)};
            data.reference(reference, false);
        }
        std::cerr << "Verified reference: degree=" << settings.degree << " fixture=" << data.name
                  << (settings.degree <= 67 ? " schoolbook/native-division\n" : " mixed-radix/identities\n");
        benchmark_case<native_backend> native(data, "mixed_radix", order);
        benchmark_case<small_backend> small(data, "small_prime", std::bit_ceil(2 * settings.degree - 1));
        for (std::size_t op = 0; op < operation_names.size(); ++op) {
            for (std::size_t trial = 0; trial < settings.warmups + settings.trials; ++trial) {
                const auto phase = trial < settings.warmups ? "warmup" : "trial";
                // Alternate the backend order to reduce systematic timing bias. Operations remain serial.
                if ((op + trial) % 2 == 0) {
                    native.sample(static_cast<operation>(op), phase, trial);
                    small.sample(static_cast<operation>(op), phase, trial);
                } else {
                    small.sample(static_cast<operation>(op), phase, trial);
                    native.sample(static_cast<operation>(op), phase, trial);
                }
            }
        }
        std::cerr << "Verified all samples: degree=" << settings.degree << " fixture=" << data.name << '\n';
    }
}    // namespace

int main(int argc, char **argv) {
    try {
        options settings;
        for (int i = 1; i < argc; ++i) {
            const std::string_view argument(argv[i]);
            if (argument == "--check") {
                settings.check = true;
            } else if (argument == "--input-dir") {
                require(++i < argc, "missing input directory");
                settings.input_directory = argv[i];
            } else {
                auto *target = argument == "--degree"  ? &settings.degree :
                               argument == "--trials"  ? &settings.trials :
                               argument == "--warmups" ? &settings.warmups :
                                                         nullptr;
                require(target && ++i < argc,
                        "usage: [--check] [--degree N] [--trials N] [--warmups N] [--input-dir DIR]");
                const std::string_view value(argv[i]);
                const auto parsed = std::from_chars(value.data(), value.data() + value.size(), *target);
                require(parsed.ec == std::errc {} && parsed.ptr == value.data() + value.size(), "invalid number");
            }
        }
        require(settings.degree >= 1 && settings.degree <= 17487 && settings.trials >= 1 && settings.trials <= 1000 &&
                    settings.warmups <= 1000,
                "require degree 1..17487, trials 1..1000, warmups 0..1000");
        std::cerr << "BN254 Fp12; boost=" << BOOST_LIB_VERSION << "; seed=" << seed
                  << "; serial; synthetic inputs; setup and checks excluded from operation samples.\n";
        std::cout << std::fixed << std::setprecision(6) << std::unitbuf
                  << "backend,degree,fixture,max_transform,phase,trial,operation,wall_ms,cpu_ms,prepared_transform,"
                     "prepared_bytes\n";
        if (settings.check) {
            settings.trials = 1;
            settings.warmups = 0;
            for (std::size_t degree : {1, 2, 4, 5, 6, 8, 9, 17, 32, 67}) {
                settings.degree = degree;
                run(settings, false);
                run(settings, true);
            }
        } else {
            run(settings, false);
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
