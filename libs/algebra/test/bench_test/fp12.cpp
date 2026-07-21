//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

#include <boost/random/mersenne_twister.hpp>

#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

namespace fields = nil::crypto3::algebra::fields;

using alt_bn128_254_fp12 = fields::fp12_2over3over2<fields::alt_bn128<254>>;
using bls12_377_fp12 = fields::fp12_2over3over2<fields::bls12<377>>;
using bls12_381_fp12 = fields::fp12_2over3over2<fields::bls12<381>>;

struct bench_config {
    std::size_t iters = 100'000ULL;
    std::size_t poolN = 1024ULL;
    std::size_t warmup = 10'000ULL;
    std::size_t samples = 20ULL;
};

static std::uint64_t now_ns() {
    using clock = std::chrono::steady_clock;
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now().time_since_epoch()).count());
}

static void do_not_optimize(void const *value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(value) : "memory");
#else
    static volatile void const *sink;
    sink = value;
#endif
}

struct bench_result {
    double total_ns = 0.0;
    double ns_per = 0.0;
    double stddev_ns_per = 0.0;
    std::size_t samples = 0;
};

template<typename Body>
bench_result run_stage(std::size_t iters, std::size_t warmup, std::size_t samples, Body &&body) {
    for (std::size_t i = 0; i < warmup; ++i) {
        body(i);
    }

    const std::size_t sample_count = (iters < samples) ? iters : samples;
    const std::size_t base_iters = iters / sample_count;
    const std::size_t extra_iters = iters % sample_count;
    std::size_t iter = 0;
    double mean_ns_per = 0.0;
    double m2_ns_per = 0.0;

    bench_result result;
    result.samples = sample_count;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const std::size_t sample_iters = base_iters + (sample < extra_iters ? 1 : 0);
        const std::uint64_t t0 = now_ns();
        for (std::size_t i = 0; i < sample_iters; ++i) {
            body(iter++);
        }
        const std::uint64_t t1 = now_ns();

        const double sample_total_ns = static_cast<double>(t1 - t0);
        const double sample_ns_per = sample_total_ns / static_cast<double>(sample_iters);
        result.total_ns += sample_total_ns;

        const double delta = sample_ns_per - mean_ns_per;
        mean_ns_per += delta / static_cast<double>(sample + 1);
        const double delta2 = sample_ns_per - mean_ns_per;
        m2_ns_per += delta * delta2;
    }

    result.ns_per = result.total_ns / static_cast<double>(iters);
    result.stddev_ns_per = sample_count > 1 ? std::sqrt(m2_ns_per / static_cast<double>(sample_count - 1)) : 0.0;
    return result;
}

static void print_stage(std::string_view operation, const bench_result &result) {
    std::cout << "  " << std::left << std::setw(12) << operation << " per=" << std::setw(12) << std::fixed
              << std::setprecision(2) << result.ns_per << " ns"
              << " stddev=" << std::setw(12) << result.stddev_ns_per << " ns\n";
}

template<typename Fp12Field>
static void run_field_benchmark(const std::string &name, const bench_config &config) {
    using base_field_type = typename Fp12Field::base_field_type;
    using value_type = typename Fp12Field::value_type;

    boost::random::mt19937 rng(1234);

    std::vector<value_type> xs(config.poolN);
    std::vector<value_type> ys(config.poolN);

    for (std::size_t i = 0; i < config.poolN; ++i) {
        xs[i] = nil::crypto3::algebra::random_element<Fp12Field>(rng);
        ys[i] = nil::crypto3::algebra::random_element<Fp12Field>(rng);
    }

    value_type acc;

    std::cout << "\n" << name << " fp12 benchmark\n";
    std::cout << "  modulus_bits=" << base_field_type::modulus_bits << " arity=" << Fp12Field::arity << "\n";

    print_stage("Fp12 mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    acc = xs[idx] * ys[idx];
                    do_not_optimize(&acc);
                }));
}

static void print_usage(const char *argv0) {
    std::cerr << "usage: " << argv0 << " [iters>0] [poolN>0] [warmup] [samples>0]\n";
}

int main(int argc, char **argv) {
    bench_config config;
    config.iters = (argc >= 2) ? std::strtoull(argv[1], nullptr, 10) : config.iters;
    config.poolN = (argc >= 3) ? std::strtoull(argv[2], nullptr, 10) : config.poolN;
    config.warmup = (argc >= 4) ? std::strtoull(argv[3], nullptr, 10) : config.warmup;
    config.samples = (argc >= 5) ? std::strtoull(argv[4], nullptr, 10) : config.samples;

    if (config.iters == 0 || config.poolN == 0 || config.samples == 0) {
        print_usage(argv[0]);
        return 1;
    }

    std::cout << "Fp12 multiplication benchmark (crypto3)\n";
    std::cout << "iters=" << config.iters << " poolN=" << config.poolN << " warmup=" << config.warmup
              << " samples=" << config.samples << "\n";

    run_field_benchmark<alt_bn128_254_fp12>("alt_bn128_254_fp12", config);
    run_field_benchmark<bls12_377_fp12>("bls12_377_fp12", config);
    run_field_benchmark<bls12_381_fp12>("bls12_381_fp12", config);

    return 0;
}
