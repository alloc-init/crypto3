//---------------------------------------------------------------------------//
// Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
//
// SPDX-License-Identifier: MIT

#include <iostream>
#include <vector>
#include <chrono>

#include <nil/crypto3/algebra/random_element.hpp>

#define CRYPTO3_RUN_BENCHMARK(bench_name, results_name, bench_type_A, bench_type_B, bench_type_C, expression) \
    do { \
        std::cout << "# " << bench_name << " " /*<< std::endl*/; \
        using duration = std::chrono::duration<double, std::nano>; \
        std::size_t long_batch; \
        { \
            bench_type_A::type_name::value_type A; bench_type_B::type_name::value_type B; bench_type_C::type_name::value_type C; \
            B = bench_type_B::sample_random_element(); \
            C = bench_type_C::sample_random_element(); \
            auto start = std::chrono::high_resolution_clock::now(); \
            std::size_t ESTIMATION_BATCH = 1000; \
            for(std::size_t i = 0; i < ESTIMATION_BATCH; ++i) { expression; } \
            auto finish = std::chrono::high_resolution_clock::now(); \
            duration sample = finish - start; \
            std::size_t WARMUP = 5; \
            std::size_t warmup_batch = 1e9*ESTIMATION_BATCH*WARMUP/sample.count(); \
            /*std::cout << "Warming up.." << std::endl;*/ \
            for(std::size_t i = 0; i < warmup_batch; ++i) { expression; } \
            long_batch = 5*1e9*ESTIMATION_BATCH/ESTIMATION_BATCH/sample.count(); \
        } \
        /* std::cout << "long_batch: " << long_batch << std::endl; */ \
        /* std::cout << "Measuring.. " << std::endl; */ \
        constexpr std::size_t MEASUREMENTS = 100; \
        constexpr std::size_t START_MULT = 2, FINISH_MULT = 10; \
        std::array<std::pair<std::size_t, duration>, (FINISH_MULT - START_MULT) * MEASUREMENTS> durs; \
        { \
            /* CAPACITY_X - create arrays that do not fit into L1 cache */ \
            constexpr std::size_t CAPACITY_B = 32768 /*sysconf(_SC_LEVEL1_DCACHE_SIZE)*/*2/sizeof(bench_type_B::type_name::value_type); \
            constexpr std::size_t CAPACITY_C = 32768 /*sysconf(_SC_LEVEL1_DCACHE_SIZE)*/*2/sizeof(bench_type_C::type_name::value_type); \
            /* STRIDE_X - stride through array so each time it misses the cache line */ \
            std::size_t STRIDE_B = 1 + sysconf(_SC_LEVEL1_DCACHE_LINESIZE)*2 / sizeof(bench_type_B::type_name::value_type); \
            std::size_t STRIDE_C = 1 + sysconf(_SC_LEVEL1_DCACHE_LINESIZE)*2 / sizeof(bench_type_C::type_name::value_type); \
            std::array<bench_type_B::type_name::value_type, CAPACITY_B> B_array; \
            std::array<bench_type_C::type_name::value_type, CAPACITY_C> C_array; \
            for(std::size_t i = 0 ; i < CAPACITY_B; ++i ) { B_array[i] = bench_type_B::sample_random_element(); } \
            for(std::size_t i = 0 ; i < CAPACITY_C; ++i ) { C_array[i] = bench_type_C::sample_random_element(); } \
            /*std::cout << "Random elements sampled.." << std::endl;*/ \
            std::vector<bench_type_A::type_name::value_type> acc(long_batch); \
            std::size_t B_idx = 0, C_idx = 0; \
            std::size_t m = 0; \
            for(std::size_t d = START_MULT; d < FINISH_MULT; ++d) { \
                for(std::size_t b = 0 ; b < MEASUREMENTS; ++b) { \
                    auto start = std::chrono::high_resolution_clock::now(); \
                    for(std::size_t i = 0; i < d*long_batch; ++i) { \
                        bench_type_A::type_name::value_type &A = acc[i % long_batch]; \
                        bench_type_B::type_name::value_type &B = B_array[B_idx % CAPACITY_B]; \
                        bench_type_C::type_name::value_type &C = C_array[C_idx % CAPACITY_C]; \
                        expression; \
                        B_idx += STRIDE_B; C_idx += STRIDE_C; \
                    } \
                    auto finish = std::chrono::high_resolution_clock::now(); \
                    durs[m++] = std::make_pair(d*long_batch, finish - start); \
                } \
                /*std::cout << "batch " << d << ":" << std::fixed << std::setprecision(3) << durs[d].count() << std::endl;*/ \
            } \
        } \
        double xy = 0, x2 = 0; \
        for(auto const& p: durs) { \
            xy += p.first*p.second.count(); x2 += p.first*p.first; \
        } \
        double b = xy/x2; \
        std::cout << std::fixed << std::setprecision(3) << b << std::endl; \
    } while(0);

#if 0
std::sort(durs.begin(), durs.end()); \
        double avg = 0, stdiv = 0, m2 = 0; \
        for(std::size_t i = MEASUREMENTS*0.1; i < MEASUREMENTS*0.9; ++i ) { \
            avg += durs[i].count(); \
            m2 += durs[i].count()*durs[i].count(); \
        } \
        avg /= MEASUREMENTS*0.8; \
        stdiv = sqrt(m2/(MEASUREMENTS*0.8) - avg*avg); \
        std::cout << " median: "  << std::fixed << std::setprecision(2) << durs[MEASUREMENTS/2].count() << "ns "; \
        std::cout << "mean: "  << std::fixed << std::setprecision(2) << avg << "ns "; \
        std::cout << "stdiv: " << std::fixed << std::setprecision(2) << stdiv << std::endl; \

#endif

template<typename A>
class bench_type {
    public:
        typedef A type_name;
        static typename A::value_type sample_random_element() {
            return nil::crypto3::algebra::random_element<A>();
        }
};
