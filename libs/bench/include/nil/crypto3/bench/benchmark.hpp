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
        std::cout << "Running test: " << bench_name << std::endl; \
        std::cout << "Estimating batch sizes.." << std::endl; \
        { \
            bench_type_A::type_name::value_type A; bench_type_B::type_name::value_type B; bench_type_C::type_name::value_type C; \
            B = bench_type_B::sample_random_element(); \
            C = bench_type_C::sample_random_element(); \
            using duration = std::chrono::duration<double, std::nano>; \
            auto start = std::chrono::high_resolution_clock::now(); \
            for(int i = 0; i < 10; ++i) { expression; } \
            auto finish = std::chrono::high_resolution_clock::now(); \
            duration sample = finish - start; \
            std::cout << "batch of 10 runs: " << sample.count() << "ns" << std::endl; \
        } \
    } while(0);

template<typename A>
class bench_type {
    public:
        typedef A type_name;
        static typename A::value_type sample_random_element() {
            return nil::crypto3::algebra::random_element<A>();
        }
};
