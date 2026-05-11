//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

using curve_type = nil::crypto3::algebra::curves::alt_bn128<254>;
using fp12_type = typename curve_type::gt_type;
using fp12_value_type = typename fp12_type::value_type;

static std::uint64_t now_ns() {
    using clock = std::chrono::steady_clock;
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            clock::now().time_since_epoch())
            .count());
}

static void do_not_optimize(void const* value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(value) : "memory");
#else
    static volatile void const* sink;
    sink = value;
#endif
}

int main(int argc, char** argv) {
    std::size_t iters = (argc >= 2) ? std::strtoull(argv[1], nullptr, 10) : 100'000ULL;
    std::size_t poolN = (argc >= 3) ? std::strtoull(argv[2], nullptr, 10) : 1024ULL;
    std::size_t warmup = (argc >= 4) ? std::strtoull(argv[3], nullptr, 10) : 1'000'000ULL;

    if (iters == 0 || poolN == 0) {
        std::cerr << "usage: " << argv[0] << " [iters>0] [poolN>0] [warmup]\n";
        return 1;
    }

    std::vector<fp12_value_type> xs(poolN);
    std::vector<fp12_value_type> ys(poolN);

    for (std::size_t i = 0; i < poolN; ++i) {
        xs[i] = nil::crypto3::algebra::random_element<fp12_type>();
        ys[i] = nil::crypto3::algebra::random_element<fp12_type>();
    }

    fp12_value_type acc;
    for (std::size_t i = 0; i < warmup; ++i) {
        const fp12_value_type& a = xs[i % poolN];
        const fp12_value_type& b = ys[i % poolN];
        acc = a * b;
        do_not_optimize(&acc);
    }

    const std::uint64_t t0 = now_ns();
    for (std::size_t i = 0; i < iters; ++i) {
        const fp12_value_type& a = xs[i % poolN];
        const fp12_value_type& b = ys[i % poolN];
        acc = a * b;
        do_not_optimize(&acc);
    }
    const std::uint64_t t1 = now_ns();

    const double total_ns = static_cast<double>(t1 - t0);
    const double ns_per = total_ns / static_cast<double>(iters);

    std::cout << "BN254 Fp12 mul benchmark (crypto3)\n";
    std::cout << "iters=" << iters << " poolN=" << poolN << " warmup=" << warmup << "\n";
    std::cout << "total: " << total_ns * 1e-9 << " s\n";
    std::cout << "per mul: " << ns_per << " ns\n";
    std::cout << "throughput: " << (1e9 / ns_per) << " mul/s\n";
    std::cout << "acc: " << acc << "\n";

    return 0;
}
