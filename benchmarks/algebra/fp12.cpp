//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#include <chrono>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <boost/random/mersenne_twister.hpp>

#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

namespace fields = nil::crypto3::algebra::fields;
namespace fp12_limb_ops = nil::crypto3::algebra::fields::detail::alt_bn128_fp12_limb_ops;

using alt_bn128_254_fp12 = fields::fp12_2over3over2<fields::alt_bn128<254>>;
using bn128_254_fp12 = fields::fp12_2over3over2<fields::bn128<254>>;
using bls12_377_fp12 = fields::fp12_2over3over2<fields::bls12<377>>;
using bls12_381_fp12 = fields::fp12_2over3over2<fields::bls12<381>>;

template<typename FieldType>
struct fp12_field_descriptor;

template<>
struct fp12_field_descriptor<alt_bn128_254_fp12> {
    constexpr static std::string_view name = "alt_bn128_254";
    constexpr static std::uint32_t seed = 0x25401;
};

template<>
struct fp12_field_descriptor<bn128_254_fp12> {
    constexpr static std::string_view name = "bn128_254";
    constexpr static std::uint32_t seed = 0x25402;
};

template<>
struct fp12_field_descriptor<bls12_377_fp12> {
    constexpr static std::string_view name = "bls12_377";
    constexpr static std::uint32_t seed = 0x37712;
};

template<>
struct fp12_field_descriptor<bls12_381_fp12> {
    constexpr static std::string_view name = "bls12_381";
    constexpr static std::uint32_t seed = 0x38112;
};

struct bench_config {
    std::size_t iters = 100'000ULL;
    std::size_t poolN = 1024ULL;
    std::size_t warmup = 10'000ULL;
    std::size_t samples = 20ULL;
};

template<typename Fp12Field>
struct has_four_limb_backend {
    using backend_type = typename Fp12Field::base_field_type::modular_backend;

    constexpr static bool value = backend_type::limb_bits == fp12_limb_ops::limb_bits &&
                                  backend_type::internal_limb_count == fp12_limb_ops::base_value_limb_count;
};

template<typename Fp12Field>
concept has_limb_fp12_fast_stages =
    has_four_limb_backend<Fp12Field>::value &&
    requires(typename Fp12Field::value_type x, typename Fp12Field::value_type y) {
        typename Fp12Field::extension_policy::base_value_type;
        { Fp12Field::extension_policy::multiply(x, y) } -> std::same_as<typename Fp12Field::value_type>;
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

struct stage_result {
    std::string field;
    std::string name;
    bench_result result;
};

static std::vector<stage_result> &stage_results() {
    static std::vector<stage_result> results;
    return results;
}

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

static void print_stage(std::string_view field, std::string name, const bench_result &result) {
    std::cout << "  " << std::left << std::setw(28) << name << " per=" << std::setw(12) << std::fixed
              << std::setprecision(2) << result.ns_per << " ns"
              << " stddev=" << std::setw(12) << result.stddev_ns_per << " ns\n";
    stage_results().push_back({std::string(field), std::move(name), result});
}

static void print_csv_results() {
    std::cout << "\nCSV\n";
    std::cout << "Field,Operation,crypto3 (ns),stddev (ns)\n";
    for (const stage_result &stage : stage_results()) {
        std::cout << stage.field << "," << stage.name << "," << std::fixed << std::setprecision(2)
                  << stage.result.ns_per << "," << stage.result.stddev_ns_per << "\n";
    }
}

template<typename Fp12Value, typename Fp6Base>
struct fp12_packed_parts {
    Fp6Base a;
    Fp6Base b;
    Fp6Base c;
    Fp6Base d;

    fp12_packed_parts(const Fp12Value &x, const Fp12Value &y) :
        a(x.data[0]), b(x.data[1]), c(y.data[0]), d(y.data[1]) {
    }
};

template<typename Fp12Value, typename Fp6Dbl, typename PackedParts>
static void multiply_prepacked_fp12(Fp12Value &ret, const PackedParts &parts) {
    Fp6Dbl ac, bd, z;
    Fp6Dbl::mul_pre(ac, parts.a, parts.c);
    Fp6Dbl::mul_pre(bd, parts.b, parts.d);

    Fp6Dbl::mul_pre(z, parts.a + parts.b, parts.c + parts.d);
    z -= ac;
    z -= bd;
    z.to_underlying(ret.data[1]);

    Fp6Dbl::mul_v_add(z, bd, ac);
    z.to_underlying(ret.data[0]);
}

template<typename Fp12Field>
static void run_limb_fast_stages(const bench_config &config,
                                 const std::vector<typename Fp12Field::base_field_type::value_type> &fpxs,
                                 const std::vector<typename Fp12Field::base_field_type::value_type> &fpys,
                                 const std::vector<typename Fp12Field::extension_policy::non_residue_field_type::value_type>
                                     &fp2xs,
                                 const std::vector<typename Fp12Field::extension_policy::non_residue_field_type::value_type>
                                     &fp2ys,
                                 const std::vector<typename Fp12Field::extension_policy::underlying_field_type::value_type>
                                     &fp6xs,
                                 const std::vector<typename Fp12Field::extension_policy::underlying_field_type::value_type>
                                     &fp6ys,
                                 const std::vector<typename Fp12Field::value_type> &xs,
                                 const std::vector<typename Fp12Field::value_type> &ys)
    requires has_limb_fp12_fast_stages<Fp12Field>
{
    constexpr std::string_view field = fp12_field_descriptor<Fp12Field>::name;

    using base_field_type = typename Fp12Field::base_field_type;
    using base_value_type = typename base_field_type::value_type;
    using fp2_value_type = typename Fp12Field::extension_policy::non_residue_field_type::value_type;
    using fp6_value_type = typename Fp12Field::extension_policy::underlying_field_type::value_type;
    using fp12_value_type = typename Fp12Field::value_type;
    using fp12_policy_type = typename Fp12Field::extension_policy;
    using fp12_fast_type =
        nil::crypto3::algebra::fields::detail::alt_bn128_fp12_fast_multiply<base_field_type, fp12_policy_type>;
    using fp2_base_type = typename fp12_fast_type::fp2_base;
    using fp2_dbl_type = typename fp12_fast_type::fp2_dbl;
    using fp6_base_type = typename fp12_fast_type::fp6_base;
    using fp6_dbl_type = typename fp12_fast_type::fp6_dbl;
    using packed_parts_type = fp12_packed_parts<fp12_value_type, fp6_base_type>;

    std::vector<typename fp12_fast_type::limb_array> fp_limbs_x(config.poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_limbs_y(config.poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_products(config.poolN);
    std::vector<typename fp12_fast_type::limb_array> fp_sum_products(config.poolN);
    std::vector<fp2_base_type> fp2_base_x(config.poolN);
    std::vector<fp2_base_type> fp2_base_y(config.poolN);
    std::vector<fp2_dbl_type> fp2_dbl_x(config.poolN);
    std::vector<fp2_dbl_type> fp2_dbl_y(config.poolN);
    std::vector<fp6_base_type> fp6_base_x(config.poolN);
    std::vector<fp6_base_type> fp6_base_y(config.poolN);
    std::vector<fp6_dbl_type> fp6_dbl_x(config.poolN);
    std::vector<fp6_dbl_type> fp6_dbl_y(config.poolN);
    std::vector<packed_parts_type> fp12_packed_parts_x;
    fp12_packed_parts_x.reserve(config.poolN);

    for (std::size_t i = 0; i < config.poolN; ++i) {
        fp_limbs_x[i] = fp12_limb_ops::load_limbs(fpxs[i].data.backend().base_data());
        fp_limbs_y[i] = fp12_limb_ops::load_limbs(fpys[i].data.backend().base_data());
        fp2_base_x[i] = fp2_base_type(fp2xs[i]);
        fp2_base_y[i] = fp2_base_type(fp2ys[i]);
        fp6_base_x[i] = fp6_base_type(fp6xs[i]);
        fp6_base_y[i] = fp6_base_type(fp6ys[i]);
        fp12_packed_parts_x.emplace_back(xs[i], ys[i]);
    }

    for (std::size_t i = 0; i < config.poolN; ++i) {
        const std::size_t next = (i + 1) % config.poolN;
        fp12_limb_ops::multiply_4x4(fp_products[i].data(), fp_limbs_x[i].data(), fp_limbs_y[i].data());
        fp12_limb_ops::multiply_4x4(fp_sum_products[i].data(), fp_limbs_y[i].data(), fp_limbs_x[next].data());

        fp2_dbl_type::mul_pre(fp2_dbl_x[i], fp2_base_x[i], fp2_base_y[i]);
        fp2_dbl_type::mul_pre(fp2_dbl_y[i], fp2_base_y[i], fp2_base_x[next]);

        fp6_dbl_type::mul_pre(fp6_dbl_x[i], fp6_base_x[i], fp6_base_y[i]);
        fp6_dbl_type::mul_pre(fp6_dbl_y[i], fp6_base_y[i], fp6_base_x[next]);
    }

    base_value_type fp_acc;
    typename fp12_fast_type::limb_array fp_limb_acc;
    typename fp12_fast_type::limb_array fp_dbl_acc;
    fp2_dbl_type fp2_pre_acc;
    fp2_value_type fp2_acc;
    fp6_dbl_type fp6_pre_acc;
    fp6_value_type fp6_acc;
    fp12_value_type fp12_acc;

    print_stage(field, "Fast Fp limb 4x4", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp12_limb_ops::multiply_4x4(fp_limb_acc.data(), fp_limbs_x[idx].data(), fp_limbs_y[idx].data());
                    do_not_optimize(&fp_limb_acc);
                }));

    print_stage(field, "Fast Fp limb REDC", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp12_limb_ops::montgomery_reduce<base_field_type>(fp_limb_acc.data(), fp_products[idx]);
                    do_not_optimize(&fp_limb_acc);
                }));

    print_stage(field, "Fast Fp dbl reduce", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp12_limb_ops::montgomery_reduce<base_field_type>(
                        reinterpret_cast<typename fp12_fast_type::limb *>(fp_acc.data.backend().base_data().limbs()),
                        fp_products[idx]);
                    do_not_optimize(&fp_acc);
                }));

    print_stage(field, "Fast Fp dbl copy", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp_dbl_acc = fp_products[idx];
                    do_not_optimize(&fp_dbl_acc);
                }));

    print_stage(field, "Fast Fp dbl add", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp12_limb_ops::add_8_limbs_mod<base_field_type>(fp_dbl_acc, fp_products[idx],
                                                                    fp_sum_products[idx]);
                    do_not_optimize(&fp_dbl_acc);
                }));

    print_stage(field, "Fast Fp dbl sub", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp12_limb_ops::subtract_8_limbs_mod<base_field_type>(fp_dbl_acc, fp_products[idx],
                                                                         fp_sum_products[idx]);
                    do_not_optimize(&fp_dbl_acc);
                }));

    print_stage(field, "Fast Fp dbl mul_by_9",
                run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp12_limb_ops::mul_8_limbs_by_9<base_field_type>(fp_dbl_acc, fp_products[idx]);
                    do_not_optimize(&fp_dbl_acc);
                }));

    print_stage(field, "Fast Fp2 pack", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    const fp2_base_type x(fp2xs[idx]);
                    do_not_optimize(&x);
                }));

    print_stage(field, "Fast Fp2 pre mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp2_dbl_type::mul_pre(fp2_pre_acc, fp2_base_x[idx], fp2_base_y[idx]);
                    do_not_optimize(&fp2_pre_acc);
                }));

    print_stage(field, "Fast Fp2 reduce", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp2_dbl_x[idx].to_non_residue(fp2_acc);
                    do_not_optimize(&fp2_acc);
                }));

    print_stage(field, "Fast Fp2 lazy mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp2_dbl_type::mul_pre(fp2_pre_acc, fp2_base_x[idx], fp2_base_y[idx]);
                    fp2_pre_acc.to_non_residue(fp2_acc);
                    do_not_optimize(&fp2_acc);
                }));

    print_stage(field, "Fast Fp2 dbl add", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp2_pre_acc = fp2_dbl_x[idx];
                    fp2_pre_acc += fp2_dbl_y[idx];
                    do_not_optimize(&fp2_pre_acc);
                }));

    print_stage(field, "Fast Fp2 dbl sub", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp2_pre_acc = fp2_dbl_x[idx];
                    fp2_pre_acc -= fp2_dbl_y[idx];
                    do_not_optimize(&fp2_pre_acc);
                }));

    print_stage(field, "Fast Fp6 pack", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    const fp6_base_type x(fp6xs[idx]);
                    do_not_optimize(&x);
                }));

    print_stage(field, "Fast Fp6 pre mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp6_dbl_type::mul_pre(fp6_pre_acc, fp6_base_x[idx], fp6_base_y[idx]);
                    do_not_optimize(&fp6_pre_acc);
                }));

    print_stage(field, "Fast Fp6 pre sum", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    const std::size_t next = (idx + 1) % config.poolN;
                    const fp6_base_type x_sum = fp6_base_x[idx] + fp6_base_y[idx];
                    const fp6_base_type y_sum = fp6_base_x[next] + fp6_base_y[next];
                    fp6_dbl_type::mul_pre(fp6_pre_acc, x_sum, y_sum);
                    do_not_optimize(&fp6_pre_acc);
                }));

    print_stage(field, "Fast Fp6 reduce", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp6_dbl_x[idx].to_underlying(fp6_acc);
                    do_not_optimize(&fp6_acc);
                }));

    print_stage(field, "Fast Fp6 lazy mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp6_dbl_type::mul_pre(fp6_pre_acc, fp6_base_x[idx], fp6_base_y[idx]);
                    fp6_pre_acc.to_underlying(fp6_acc);
                    do_not_optimize(&fp6_acc);
                }));

    print_stage(field, "Fast Fp6 dbl add", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp6_pre_acc = fp6_dbl_x[idx];
                    fp6_pre_acc += fp6_dbl_y[idx];
                    do_not_optimize(&fp6_pre_acc);
                }));

    print_stage(field, "Fast Fp6 dbl sub", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp6_pre_acc = fp6_dbl_x[idx];
                    fp6_pre_acc -= fp6_dbl_y[idx];
                    do_not_optimize(&fp6_pre_acc);
                }));

    print_stage(field, "Fast Fp12 pack inputs",
                run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    const packed_parts_type parts(xs[idx], ys[idx]);
                    do_not_optimize(&parts);
                }));

    print_stage(field, "Fast Fp12 core prepacked",
                run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    multiply_prepacked_fp12<fp12_value_type, fp6_dbl_type>(fp12_acc, fp12_packed_parts_x[idx]);
                    do_not_optimize(&fp12_acc);
                }));
}

template<typename Fp12Field>
static void run_field_benchmark(const bench_config &config) {
    constexpr std::string_view field = fp12_field_descriptor<Fp12Field>::name;

    using base_field_type = typename Fp12Field::base_field_type;
    using fp2_type = typename Fp12Field::extension_policy::non_residue_field_type;
    using fp6_type = typename Fp12Field::extension_policy::underlying_field_type;
    using base_value_type = typename base_field_type::value_type;
    using fp2_value_type = typename fp2_type::value_type;
    using fp6_value_type = typename fp6_type::value_type;
    using fp12_value_type = typename Fp12Field::value_type;

    boost::random::mt19937 rng(fp12_field_descriptor<Fp12Field>::seed);

    std::vector<base_value_type> fpxs(config.poolN);
    std::vector<base_value_type> fpys(config.poolN);
    std::vector<fp2_value_type> fp2xs(config.poolN);
    std::vector<fp2_value_type> fp2ys(config.poolN);
    std::vector<fp6_value_type> fp6xs(config.poolN);
    std::vector<fp6_value_type> fp6ys(config.poolN);
    std::vector<fp12_value_type> xs(config.poolN);
    std::vector<fp12_value_type> ys(config.poolN);

    for (std::size_t i = 0; i < config.poolN; ++i) {
        fpxs[i] = nil::crypto3::algebra::random_element<base_field_type>(rng);
        fpys[i] = nil::crypto3::algebra::random_element<base_field_type>(rng);
        fp2xs[i] = nil::crypto3::algebra::random_element<fp2_type>(rng);
        fp2ys[i] = nil::crypto3::algebra::random_element<fp2_type>(rng);
        fp6xs[i] = nil::crypto3::algebra::random_element<fp6_type>(rng);
        fp6ys[i] = nil::crypto3::algebra::random_element<fp6_type>(rng);
        xs[i] = nil::crypto3::algebra::random_element<Fp12Field>(rng);
        ys[i] = nil::crypto3::algebra::random_element<Fp12Field>(rng);
    }

    base_value_type fp_acc;
    fp2_value_type fp2_acc;
    fp6_value_type fp6_acc;
    fp12_value_type fp12_acc;

    std::cout << "\n" << field << " fp12 benchmark\n";
    std::cout << "  modulus_bits=" << base_field_type::modulus_bits << " arity=" << Fp12Field::arity << "\n";

    print_stage(field, "Fp mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp_acc = fpxs[idx] * fpys[idx];
                    do_not_optimize(&fp_acc);
                }));

    print_stage(field, "Fp2 mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp2_acc = fp2xs[idx] * fp2ys[idx];
                    do_not_optimize(&fp2_acc);
                }));

    print_stage(field, "Fp6 mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp6_acc = fp6xs[idx] * fp6ys[idx];
                    do_not_optimize(&fp6_acc);
                }));

    print_stage(field, "Fp12 mul", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp12_acc = xs[idx] * ys[idx];
                    do_not_optimize(&fp12_acc);
                }));

    print_stage(field, "Fp12 square", run_stage(config.iters, config.warmup, config.samples, [&](std::size_t i) {
                    const std::size_t idx = i % config.poolN;
                    fp12_acc = xs[idx].squared();
                    do_not_optimize(&fp12_acc);
                }));

    if constexpr (has_limb_fp12_fast_stages<Fp12Field>) {
        run_limb_fast_stages<Fp12Field>(config, fpxs, fpys, fp2xs, fp2ys, fp6xs, fp6ys, xs, ys);
    } else {
        std::cout << "  Fast limb stages skipped for " << field << "\n";
    }
}

template<typename Fp12Field>
static bool run_field_if_selected(std::string_view filter, const bench_config &config) {
    constexpr std::string_view field = fp12_field_descriptor<Fp12Field>::name;
    if (filter == "all" || filter == field) {
        run_field_benchmark<Fp12Field>(config);
        return true;
    }
    return false;
}

static void print_usage(const char *argv0) {
    std::cerr << "usage: " << argv0 << " [iters>0] [poolN>0] [warmup] [samples>0] "
              << "[all|alt_bn128_254|bn128_254|bls12_377|bls12_381]\n";
}

int main(int argc, char **argv) {
    bench_config config;
    config.iters = (argc >= 2) ? std::strtoull(argv[1], nullptr, 10) : config.iters;
    config.poolN = (argc >= 3) ? std::strtoull(argv[2], nullptr, 10) : config.poolN;
    config.warmup = (argc >= 4) ? std::strtoull(argv[3], nullptr, 10) : config.warmup;
    config.samples = (argc >= 5) ? std::strtoull(argv[4], nullptr, 10) : config.samples;
    const std::string_view field_filter = (argc >= 6) ? std::string_view(argv[5]) : std::string_view("all");

    if (config.iters == 0 || config.poolN == 0 || config.samples == 0) {
        print_usage(argv[0]);
        return 1;
    }

    std::cout << "Fp12 tower mul benchmark (crypto3)\n";
    std::cout << "iters=" << config.iters << " poolN=" << config.poolN << " warmup=" << config.warmup
              << " samples=" << config.samples << " fields=" << field_filter << "\n";

    std::size_t fields_run = 0;
    fields_run += run_field_if_selected<alt_bn128_254_fp12>(field_filter, config) ? 1u : 0u;
    fields_run += run_field_if_selected<bn128_254_fp12>(field_filter, config) ? 1u : 0u;
    fields_run += run_field_if_selected<bls12_377_fp12>(field_filter, config) ? 1u : 0u;
    fields_run += run_field_if_selected<bls12_381_fp12>(field_filter, config) ? 1u : 0u;

    if (fields_run == 0) {
        std::cerr << "unknown field filter: " << field_filter << "\n";
        print_usage(argv[0]);
        return 1;
    }

    print_csv_results();

    return 0;
}
