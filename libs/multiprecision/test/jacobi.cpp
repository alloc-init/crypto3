///////////////////////////////////////////////////////////////
//  Copyright (c) 2020 Mikhail Komarov.
//  Copyright (c) 2021 Aleksei Moskvin <alalmoskvin@gmail.com>
//  Distributed under the Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

#define BOOST_TEST_MODULE jacobi_multiprecision_test

#ifdef _MSC_VER
#define _SCL_SECURE_NO_WARNINGS
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <cstddef>

#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>
#include <nil/crypto3/multiprecision/cpp_int_modular/literals.hpp>

#include <nil/crypto3/multiprecision/jacobi.hpp>

template<typename T, std::size_t N>
T decimal_number(const char (&decimal)[N]) {
    T result = 0u;
    for (std::size_t i = 0; i + 1 < N; ++i) {
        result *= 10u;
        result += static_cast<unsigned>(decimal[i] - '0');
    }
    return result;
}

template<typename T>
void test() {
    using namespace boost::multiprecision;

    BOOST_CHECK_EQUAL(jacobi(T(5u), T(9u)), 1);
    BOOST_CHECK_EQUAL(jacobi(T(1u), T(27u)), 1);
    BOOST_CHECK_EQUAL(jacobi(T(2u), T(27u)), -1);
    BOOST_CHECK_EQUAL(jacobi(T(3u), T(27u)), 0);
    BOOST_CHECK_EQUAL(jacobi(T(4u), T(27u)), 1);
    BOOST_CHECK_EQUAL(jacobi(T(506u), T(1103u)), -1);

    // new tests from algebra:
    BOOST_CHECK_EQUAL(
        jacobi(T(76749407),
               decimal_number<T>("21888242871839275222246405745257275088696311157297823662689037894645226208583")),
        -1);
    BOOST_CHECK_EQUAL(
        jacobi(T(76749407),
               decimal_number<T>("52435875175126190479447740508185965837690552500527637822603658699938581184513")),
        -1);
    BOOST_CHECK_EQUAL(
        jacobi(T(76749407),
               decimal_number<T>(
                   "184014718449470976641732519409003087090464839378677150738808371108642394980708029691298675282643534"
                   "004240329819623250370915624553421958933568063851025402776443788222357196988103580408517486311017895"
                   "169440640343141708939276039764731720833213215559801639066799283898191098079351209268491644339667178"
                   "604494222971572788971054374438281331602764950963417101448891412424011588862068850113410088177801409"
                   "279786489730635599081340855930762685458174837104230446238204727771628459008795937374640002233231333"
                   "609522446689297900090549115400764760910459967591503490110147729489296261451835450258703237412701103"
                   "14006814529932451772897")),
        -1);
}

BOOST_AUTO_TEST_SUITE(jacobi_tests)

BOOST_AUTO_TEST_CASE(jacobi_test) {
    using namespace boost::multiprecision;

    test<cpp_int>();

    constexpr auto a = 0x4931a5f_cppui_modular256;
    constexpr auto b = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001_cppui_modular256;
    static_assert(jacobi(a, b) == -1, "jacobi error");
}

BOOST_AUTO_TEST_SUITE_END()
