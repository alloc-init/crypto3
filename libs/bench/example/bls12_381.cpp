//---------------------------------------------------------------------------//
// Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
//
// SPDX-License-Identifier: MIT

#include <iostream>

#include <nil/crypto3/bench/benchmark.hpp>

#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/curves/mnt4.hpp>
#include <nil/crypto3/algebra/curves/pallas.hpp>

#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

#include <nil/crypto3/algebra/random_element.hpp>

int main()
{
    {
        using bench_type_A = bench_type<nil::crypto3::algebra::curves::mnt4_298::base_field_type>;
        using bench_type_B = bench_type_A;
        using bench_type_C = bench_type_A;

        CRYPTO3_RUN_BENCHMARK("MNT4-298 base field multiplication", "sample.json",
                bench_type_A, bench_type_B, bench_type_C,
                A = B * C
                );
    }
    {
        using bench_type_A = bench_type<nil::crypto3::algebra::curves::mnt4_298::base_field_type>;
        using bench_type_B = bench_type_A;
        using bench_type_C = bench_type_A;
        auto op = [](
                bench_type_A::type_name::value_type & la,
                bench_type_B::type_name::value_type const& lb,
                bench_type_C::type_name::value_type const& lc) { 
            la = lb * lc;
        };

        CRYPTO3_RUN_BENCHMARK("MNT4-298 base field multiplication, lambda", "sample.json",
                bench_type_A, bench_type_B, bench_type_C,

                op(A,B,C);
                );
    }
    /*
    {
        using bench_type_A = bench_type<nil::crypto3::algebra::curves::mnt4_298::g1_type<>>;
        using bench_type_B = bench_type_A;
        using bench_type_C = bench_type<nil::crypto3::algebra::curves::mnt4_298::scalar_field_type>;

        CRYPTO3_RUN_BENCHMARK("mnt4-298 Curve point addition", "sample.json",
                bench_type_A, bench_type_B, bench_type_B,
                A = B + C
                );
        CRYPTO3_RUN_BENCHMARK("mnt4-298 Curve point addition inplace", "sample.json",
                bench_type_A, bench_type_B, bench_type_B,
                A += B
                );
        CRYPTO3_RUN_BENCHMARK("mnt4-298 Curve point doubling", "sample.json",
                bench_type_A, bench_type_B, bench_type_B,
                B.double_inplace()
                );
        CRYPTO3_RUN_BENCHMARK("mnt4-298 Curve scalar multiplication", "sample.json",
                bench_type_A, bench_type_B, bench_type_C,
                A = B*C
                );
    }
    */
    /*
    {
        using bench_type_A = bench_type<nil::crypto3::algebra::curves::bls12_381::g1_type<>>;
        using bench_type_B = bench_type_A;
        using bench_type_C = bench_type<nil::crypto3::algebra::curves::bls12_381::scalar_field_type>;

        CRYPTO3_RUN_BENCHMARK("BLS12-381 Curve point addition", "sample.json",
                bench_type_A, bench_type_B, bench_type_B,
                A = B + C
                );
        CRYPTO3_RUN_BENCHMARK("BLS12-381 Curve point addition inplace", "sample.json",
                bench_type_A, bench_type_B, bench_type_B,
                A += B
                );
        CRYPTO3_RUN_BENCHMARK("BLS12-381 Curve point doubling", "sample.json",
                bench_type_A, bench_type_B, bench_type_B,
                B.double_inplace()
                );
        CRYPTO3_RUN_BENCHMARK("BLS12-381 Curve scalar multiplication", "sample.json",
                bench_type_A, bench_type_B, bench_type_C,
                A = B*C
                );
    }
    */
    /*
    {
        using bench_type_A = bench_type<nil::crypto3::algebra::curves::bls12_381::base_field_type>;
        using bench_type_B = bench_type_A;
        using bench_type_C = bench_type_A;

        CRYPTO3_RUN_BENCHMARK("BLS12-381 base field multiplication", "sample.json",
                bench_type_A, bench_type_B, bench_type_C,
                A = B * C
                );
    }
    */
    /*
    {
        using bench_type_A = bench_type<nil::crypto3::algebra::curves::bls12_381::scalar_field_type>;
        using bench_type_B = bench_type_A;
        using bench_type_C = bench_type_A;

        CRYPTO3_RUN_BENCHMARK("BLS12-381 scalar field multiplication", "sample.json",
                bench_type_A, bench_type_B, bench_type_C,
                A = B * C
                );
    }
    */
    /*
    {
        using bench_type_A = bench_type<nil::crypto3::algebra::curves::pallas::base_field_type>;
        using bench_type_B = bench_type_A;
        using bench_type_C = bench_type_A;

        CRYPTO3_RUN_BENCHMARK("Pallas base field multiplication", "sample.json",
                bench_type_A, bench_type_B, bench_type_C,
                A = B * C
                );
    }
    */
    return 0;
}
