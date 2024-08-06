//---------------------------------------------------------------------------//
// Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
//
// SPDX-License-Identifier: MIT

#ifndef CRYPTO3_BENCHMARK_HPP
#define CRYPTO3_BENCHMARK_HPP

#include <iostream>
#include <vector>
#include <array>
#include <chrono>

#include <nil/crypto3/algebra/random_element.hpp>

#define CRYPTO3_RUN_BENCHMARK(bench_name, bench_type_A, bench_type_B, bench_type_C, expression) \
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
            long_batch = 1e9*ESTIMATION_BATCH/ESTIMATION_BATCH/sample.count(); \
        } \
        /* std::cout << "long_batch: " << long_batch << std::endl; */ \
        /* std::cout << "Measuring.. " << std::endl; */ \
        constexpr std::size_t MEASUREMENTS = 100; \
        constexpr std::size_t START_MULT = 5, FINISH_MULT = 15; \
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
                        /* ### */ expression; /* ### */ \
                        B_idx += STRIDE_B; C_idx += STRIDE_C; \
                    } \
                    auto finish = std::chrono::high_resolution_clock::now(); \
                    durs[m++] = std::make_pair(d*long_batch, finish - start); \
                } \
                /*std::cout << "batch " << d << ":" << std::fixed << std::setprecision(3) << durs[d].count() << std::endl;*/ \
            } \
        } \
        std::sort(durs.begin(), durs.end(), [](auto const& a, auto const& b) { \
                return a.first == b.first ? a.second < b.second : a.first < b.first; \
                }); \
        double xy = 0, x2 = 0; \
        for(std::size_t i = 0; i < FINISH_MULT - START_MULT; ++i) { \
            for(std::size_t m = 0; m < MEASUREMENTS; ++m) { \
                /* discard top 30% outliers */ \
                if ( m < MEASUREMENTS*0.7 ) { \
                    xy += durs[i*MEASUREMENTS+m].first * durs[i*MEASUREMENTS+m].second.count(); \
                    x2 += durs[i*MEASUREMENTS+m].first * durs[i*MEASUREMENTS+m].first; \
                } \
            } \
        } \
        double b = xy/x2; \
        std::cout << std::fixed << std::setprecision(3) << b << " ns" << std::endl; \
    } while(0);

namespace nil {
    namespace crypto3 {
        namespace bench {

            template<typename A>
                class bench_type {
                    public:
                        typedef A type_name;
                        static typename A::value_type sample_random_element() {
                            return nil::crypto3::algebra::random_element<A>();
                        }
                };

        }
    }
}

#endif /* CRYPTO3_BENCHMARK_HPP */
