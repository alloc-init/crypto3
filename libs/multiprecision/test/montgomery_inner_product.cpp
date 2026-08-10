//---------------------------------------------------------------------------//
// Copyright (c) 2026 Riccardo Abbate
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#define BOOST_TEST_MODULE montgomery_inner_product_test

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <string>
#include <vector>

#include <boost/multiprecision/number.hpp>

#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>
#include <nil/crypto3/multiprecision/modular/modular_params_fixed.hpp>

using boost::multiprecision::backends::cpp_int_modular_backend;
using boost::multiprecision::backends::modular_params;

namespace {
    template<unsigned Bits>
    void check_inner_product(const char *modulus_text, std::size_t count) {
        using backend_type = cpp_int_modular_backend<Bits>;
        using number_type = boost::multiprecision::number<backend_type>;

        const number_type modulus(modulus_text);
        modular_params<backend_type> params(modulus.backend());
        std::vector<backend_type> left(count);
        std::vector<backend_type> right(count);

        backend_type expected;
        for (std::size_t i = 0; i < count; ++i) {
            // Values near the modulus maximize the raw products and exercise the accumulator's carry limb.
            left[i] = number_type(modulus - static_cast<std::size_t>(i % 17 + 1)).backend();
            right[i] = number_type(modulus - static_cast<std::size_t>(i % 23 + 1)).backend();

            backend_type product(left[i]);
            params.mod_mul(product, right[i]);
            params.mod_add(expected, product);
        }

        backend_type actual;
        params.get_mod_obj().montgomery_inner_product(actual, left.begin(), left.end(), right.begin());
        BOOST_CHECK_EQUAL(actual.compare(expected), 0);
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(montgomery_inner_product_test_suite)

BOOST_AUTO_TEST_CASE(empty_range_returns_zero) {
    check_inner_product<254>("0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47", 0);
}

BOOST_AUTO_TEST_CASE(single_product_matches_eager_multiplication) {
    check_inner_product<254>("0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47", 1);
}

BOOST_AUTO_TEST_CASE(medium_range_matches_eager_accumulation) {
    check_inner_product<254>("0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47", 64);
}

BOOST_AUTO_TEST_CASE(long_range_matches_eager_accumulation) {
    check_inner_product<255>("0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED", 1000);
}

BOOST_AUTO_TEST_CASE(full_width_modulus_matches_eager_accumulation) {
    check_inner_product<256>("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF43", 97);
}

BOOST_AUTO_TEST_CASE(trivial_backend_matches_eager_accumulation) {
    check_inner_product<31>("0x7FFFFFE7", 19);
}

BOOST_AUTO_TEST_SUITE_END()
