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
#include <string>
#include <vector>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/fields/fp2.hpp>
#include <nil/crypto3/algebra/fields/fp6_3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

using curve_type = nil::crypto3::algebra::curves::alt_bn128<254>;
using base_field_type = typename curve_type::base_field_type;
using fp2_type = nil::crypto3::algebra::fields::fp2<base_field_type>;
using fp6_type = nil::crypto3::algebra::fields::fp6_3over2<base_field_type>;
using fp12_type = typename curve_type::gt_type;

using base_value_type = typename base_field_type::value_type;
using fp2_value_type = typename fp2_type::value_type;
using fp6_value_type = typename fp6_type::value_type;
using fp12_value_type = typename fp12_type::value_type;
using fp12_policy_type = typename fp12_type::extension_policy;
using fp12_fast_type =
    nil::crypto3::algebra::fields::detail::alt_bn128_fp12_fast_multiply<base_field_type, fp12_policy_type>;
using fp2_base_type = typename fp12_fast_type::fp2_base;
using fp6_base_type = typename fp12_fast_type::fp6_base;
using fp6_dbl_type = typename fp12_fast_type::fp6_dbl;
namespace limb_ops = nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops;

static std::uint64_t now_ns() {
    using clock = std::chrono::steady_clock;
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now().time_since_epoch()).count());
}

static void do_not_optimize(void const* value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(value) : "memory");
#else
    static volatile void const* sink;
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
bench_result run_stage(std::size_t iters, std::size_t warmup, std::size_t samples, Body&& body) {
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

static void print_stage(const std::string& name, const bench_result& result) {
    std::cout << std::left << std::setw(24) << name << " per=" << std::setw(12) << result.ns_per << " ns"
              << " stddev=" << std::setw(12) << result.stddev_ns_per << " ns\n";
}

int main(int argc, char** argv) {
    std::size_t iters = (argc >= 2) ? std::strtoull(argv[1], nullptr, 10) : 100'000ULL;
    std::size_t poolN = (argc >= 3) ? std::strtoull(argv[2], nullptr, 10) : 1024ULL;
    std::size_t warmup = (argc >= 4) ? std::strtoull(argv[3], nullptr, 10) : 10'000ULL;
    std::size_t samples = (argc >= 5) ? std::strtoull(argv[4], nullptr, 10) : 20ULL;

    if (iters == 0 || poolN == 0 || samples == 0) {
        std::cerr << "usage: " << argv[0] << " [iters>0] [poolN>0] [warmup] [samples>0]\n";
        return 1;
    }

    boost::random::mt19937 rng(0);

    std::vector<base_value_type> fpxs(poolN);
    std::vector<base_value_type> fpys(poolN);
    std::vector<fp2_value_type> fp2xs(poolN);
    std::vector<fp2_value_type> fp2ys(poolN);
    std::vector<fp6_value_type> fp6xs(poolN);
    std::vector<fp6_value_type> fp6ys(poolN);
    std::vector<fp12_value_type> xs(poolN);
    std::vector<fp12_value_type> ys(poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_limbs_x(poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_limbs_y(poolN);
    std::vector<fp2_base_type> fp2_base_x(poolN);
    std::vector<fp2_base_type> fp2_base_y(poolN);
    std::vector<fp6_base_type> fp6_base_x(poolN);
    std::vector<fp6_base_type> fp6_base_y(poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_sum_x(poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_sum_y(poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_products(poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_sum_products(poolN);
    std::vector<typename fp12_fast_type::fp_dbl> fp_dbl_x(poolN);
    std::vector<typename fp12_fast_type::fp_dbl> fp_dbl_y(poolN);
    std::vector<typename fp12_fast_type::fp2_dbl> fp2_dbl_x(poolN);
    std::vector<typename fp12_fast_type::fp2_dbl> fp2_dbl_y(poolN);
    std::vector<fp6_dbl_type> fp6_dbl_x(poolN);
    std::vector<fp6_dbl_type> fp6_dbl_y(poolN);

    for (std::size_t i = 0; i < poolN; ++i) {
        fpxs[i] = nil::crypto3::algebra::random_element<base_field_type>(rng);
        fpys[i] = nil::crypto3::algebra::random_element<base_field_type>(rng);
        fp2xs[i] = nil::crypto3::algebra::random_element<fp2_type>(rng);
        fp2ys[i] = nil::crypto3::algebra::random_element<fp2_type>(rng);
        fp6xs[i] = nil::crypto3::algebra::random_element<fp6_type>(rng);
        fp6ys[i] = nil::crypto3::algebra::random_element<fp6_type>(rng);
        xs[i] = nil::crypto3::algebra::random_element<fp12_type>(rng);
        ys[i] = nil::crypto3::algebra::random_element<fp12_type>(rng);
        fp_limbs_x[i] = limb_ops::load_limbs(fpxs[i].data.backend().base_data());
        fp_limbs_y[i] = limb_ops::load_limbs(fpys[i].data.backend().base_data());
        fp2_base_x[i] = fp2_base_type(fp2xs[i]);
        fp2_base_y[i] = fp2_base_type(fp2ys[i]);
        fp6_base_x[i] = fp6_base_type(fp6xs[i]);
        fp6_base_y[i] = fp6_base_type(fp6ys[i]);
    }
    for (std::size_t i = 0; i < poolN; ++i) {
        const std::size_t next = (i + 1) % poolN;
        const std::size_t next2 = (i + 2) % poolN;
        fp_sum_x[i] = fp_limbs_x[i];
        fp_sum_y[i] = fp_limbs_y[i];
        limb_ops::add_8_limbs(fp_sum_x[i], fp_limbs_y[i]);
        limb_ops::add_8_limbs(fp_sum_x[i], fp_limbs_x[next]);
        limb_ops::add_8_limbs(fp_sum_x[i], fp_limbs_y[next]);
        limb_ops::add_8_limbs(fp_sum_y[i], fp_limbs_x[next]);
        limb_ops::add_8_limbs(fp_sum_y[i], fp_limbs_y[next]);
        limb_ops::add_8_limbs(fp_sum_y[i], fp_limbs_x[next2]);
        limb_ops::multiply_4x4(fp_products[i], fp_limbs_x[i], fp_limbs_y[i]);
        fp_dbl_x[i] = fp12_fast_type::fp_dbl(fp_products[i]);
        fp_dbl_y[i] = fp12_fast_type::fp_dbl(fp_sum_products[i]);
        fp12_fast_type::fp2_dbl::mul_pre(fp2_dbl_x[i], fp2_base_x[i], fp2_base_y[i]);
        fp12_fast_type::fp2_dbl::mul_pre(fp2_dbl_y[i], fp2_base_y[i], fp2_base_x[next]);
        fp6_dbl_type::mul_pre(fp6_dbl_x[i], fp6_base_x[i], fp6_base_y[i]);
        fp6_dbl_type::mul_pre(fp6_dbl_y[i], fp6_base_y[i], fp6_base_x[next]);
    }

    base_value_type fp_acc;
    typename fp12_fast_type::limb_array fp_limb_acc;
    typename fp12_fast_type::fp_dbl fp_dbl_acc;
    typename fp12_fast_type::fp2_dbl fp2_pre_acc;
    fp2_value_type fp2_acc;
    fp6_dbl_type fp6_pre_acc;
    fp6_value_type fp6_acc;
    fp12_value_type fp12_acc;

    std::cout << "BN254 tower mul benchmark (crypto3)\n";
    std::cout << "iters=" << iters << " poolN=" << poolN << " warmup=" << warmup << " samples=" << samples << "\n";

    print_stage("Fp limb 4x4", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    limb_ops::multiply_4x4(fp_limb_acc, fp_limbs_x[idx], fp_limbs_y[idx]);
                    do_not_optimize(&fp_limb_acc);
                }));

    print_stage("Fp dbl reduce", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    typename fp12_fast_type::fp_dbl value(fp_products[idx]);
                    value.reduce();
                    value.to_base_value(fp_acc);
                    do_not_optimize(&fp_acc);
                }));

    print_stage("Fp dbl copy", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp_dbl_acc = fp_dbl_x[idx];
                    do_not_optimize(&fp_dbl_acc);
                }));

    print_stage("Fp dbl add", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp_dbl_acc = fp_dbl_x[idx];
                    fp_dbl_acc += fp_dbl_y[idx];
                    do_not_optimize(&fp_dbl_acc);
                }));

    print_stage("Fp dbl sub", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp_dbl_acc = fp_dbl_x[idx];
                    fp_dbl_acc -= fp_dbl_y[idx];
                    do_not_optimize(&fp_dbl_acc);
                }));

    print_stage("Fp dbl mul_by_9", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp_dbl_acc = fp_dbl_x[idx];
                    fp_dbl_acc.mul_by_9();
                    do_not_optimize(&fp_dbl_acc);
                }));

    print_stage("Fp2 pre mul", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp12_fast_type::fp2_dbl fp2_pre_acc;
                    fp12_fast_type::fp2_dbl::mul_pre(fp2_pre_acc, fp2_base_x[idx], fp2_base_y[idx]);
                    do_not_optimize(&fp2_pre_acc);
                }));

    print_stage("Fp2 reduce", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    typename fp12_fast_type::fp2_dbl value(fp2_dbl_x[idx]);
                    value.reduce();
                    value.to_non_residue(fp2_acc);
                    do_not_optimize(&fp2_acc);
                }));

    print_stage("Fp2 lazy mul", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp12_fast_type::fp2_dbl fp2_pre_acc;
                    fp12_fast_type::fp2_dbl::mul_pre(fp2_pre_acc, fp2_base_x[idx], fp2_base_y[idx]);
                    fp2_pre_acc.reduce();
                    fp2_pre_acc.to_non_residue(fp2_acc);
                    do_not_optimize(&fp2_acc);
                }));

    print_stage("Fp2 dbl add", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp2_pre_acc = fp2_dbl_x[idx];
                    fp2_pre_acc += fp2_dbl_y[idx];
                    do_not_optimize(&fp2_pre_acc);
                }));

    print_stage("Fp2 dbl sub", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp2_pre_acc = fp2_dbl_x[idx];
                    fp2_pre_acc -= fp2_dbl_y[idx];
                    do_not_optimize(&fp2_pre_acc);
                }));

    print_stage("Fp6 pre mul", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp6_dbl_type::mul_pre(fp6_pre_acc, fp6_base_x[idx], fp6_base_y[idx]);
                    do_not_optimize(&fp6_pre_acc);
                }));

    print_stage("Fp6 pre sum", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    const std::size_t next = (idx + 1) % poolN;
                    const fp6_base_type x_sum = fp6_base_x[idx] + fp6_base_y[idx];
                    const fp6_base_type y_sum = fp6_base_x[next] + fp6_base_y[next];
                    fp6_dbl_type::mul_pre(fp6_pre_acc, x_sum, y_sum);
                    do_not_optimize(&fp6_pre_acc);
                }));

    print_stage("Fp6 reduce", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp6_dbl_type value(fp6_dbl_x[idx]);
                    value.reduce();
                    value.to_underlying(fp6_acc);
                    do_not_optimize(&fp6_acc);
                }));

    print_stage("Fp6 lazy mul", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp6_dbl_type::mul_pre(fp6_pre_acc, fp6_base_x[idx], fp6_base_y[idx]);
                    fp6_pre_acc.reduce();
                    fp6_pre_acc.to_underlying(fp6_acc);
                    do_not_optimize(&fp6_acc);
                }));

    print_stage("Fp6 dbl add", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp6_pre_acc = fp6_dbl_x[idx];
                    fp6_pre_acc += fp6_dbl_y[idx];
                    do_not_optimize(&fp6_pre_acc);
                }));

    print_stage("Fp6 dbl sub", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp6_pre_acc = fp6_dbl_x[idx];
                    fp6_pre_acc -= fp6_dbl_y[idx];
                    do_not_optimize(&fp6_pre_acc);
                }));

    print_stage("Fp12 mul", run_stage(iters, warmup, samples, [&](std::size_t i) {
                    const std::size_t idx = i % poolN;
                    fp12_acc = xs[idx] * ys[idx];
                    do_not_optimize(&fp12_acc);
                }));

    std::cout << "acc: " << fp12_acc << "\n";

    return 0;
}
