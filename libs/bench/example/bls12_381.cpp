//---------------------------------------------------------------------------//
// Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
//
// SPDX-License-Identifier: MIT

#include <iostream>

#include <nil/crypto3/bench/benchmark.hpp>
#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/fields/bls12/base_field.hpp>
#include <nil/crypto3/algebra/fields/bls12/scalar_field.hpp>

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/random_device.hpp>

int main()
{
    std::cout << "Hello, world" << std::endl;

    using bench_type_A = bench_type<nil::crypto3::algebra::curves::bls12_381::g1_type<>>;
    using bench_type_B = bench_type_A;
    using bench_type_C = bench_type_A; //bench_type<nil::crypto3::algebra::curves::bls12_381::scalar_field_type>;
//    using bench_type_C = bench_type<nil::crypto3::algebra::curves::bls12_381::scalar_field_type>;

    CRYPTO3_RUN_BENCHMARK("BLS12-381 Curve scalar multiplication", "sample.json",
            bench_type_A, bench_type_B, bench_type_C,
            A = B+C
            );

    return 0;
}
