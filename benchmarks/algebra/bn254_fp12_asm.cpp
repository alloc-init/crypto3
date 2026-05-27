//---------------------------------------------------------------------------//
// Copyright (c) 2026
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#include <array>
#include <cstddef>

#include <boost/random/mersenne_twister.hpp>

#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/random_element.hpp>

using curve_type = nil::crypto3::algebra::curves::alt_bn128<254>;
using fp12_type = typename curve_type::gt_type;
using fp12_value_type = typename fp12_type::value_type;

#if defined(__GNUC__) || defined(__clang__)
#define CRYPTO3_NOINLINE __attribute__((noinline))
#define CRYPTO3_USED __attribute__((used))
#else
#define CRYPTO3_NOINLINE
#define CRYPTO3_USED
#endif

static void do_not_optimize(void const *value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(value) : "memory");
#else
    static volatile void const *sink;
    sink = value;
#endif
}

extern "C" CRYPTO3_NOINLINE CRYPTO3_USED void crypto3_bn254_fp12_mul_once(
    fp12_value_type *out, const fp12_value_type *x, const fp12_value_type *y) {
    *out = (*x) * (*y);
}

extern "C" CRYPTO3_NOINLINE CRYPTO3_USED void crypto3_bn254_fp12_mul_chain(
    fp12_value_type *out, const fp12_value_type *a, const fp12_value_type *b, const fp12_value_type *c,
    const fp12_value_type *d) {
    fp12_value_type acc = (*a) * (*b);
    acc = acc * (*c);
    acc = acc * (*d);
    *out = acc;
}

int main() {
    boost::random::mt19937 rng(0);
    std::array<fp12_value_type, 4> values;
    fp12_value_type out;

    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = nil::crypto3::algebra::random_element<fp12_type>(rng);
    }

    crypto3_bn254_fp12_mul_once(&out, &values[0], &values[1]);
    do_not_optimize(&out);

    crypto3_bn254_fp12_mul_chain(&out, &out, &values[1], &values[2], &values[3]);
    do_not_optimize(&out);

    return 0;
}

#undef CRYPTO3_NOINLINE
#undef CRYPTO3_USED
