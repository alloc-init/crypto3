//---------------------------------------------------------------------------//
// MIT License
//---------------------------------------------------------------------------//

#define BOOST_TEST_MODULE matrix_test

#include <boost/test/unit_test.hpp>

#include <type_traits>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/math/matrix/backends/ublas/compressed.hpp>
#include <nil/crypto3/math/matrix/backends/ublas/regular.hpp>

namespace math = nil::crypto3::math;
namespace fields = nil::crypto3::algebra::fields;

using bn254_fp12_field = fields::fp12_2over3over2<fields::alt_bn128<254>>;
using bn254_fp12_value = typename bn254_fp12_field::value_type;

template<typename Left, typename Right>
concept LvalueAddable = requires(Left &left, Right &right) {
    left + right;
};

template<typename Left, typename Right>
concept RvalueAddable = requires {
    Left() + Right();
};

static_assert(math::MatrixBackend<math::backends::ublas::regular_matrix<int>>);
static_assert(math::MatrixBackend<math::backends::ublas::compressed_matrix<int>>);
static_assert(math::VectorBackend<math::backends::ublas::regular_vector<int>>);
static_assert(math::VectorBackend<math::backends::ublas::compressed_vector<int>>);
static_assert(math::ResizableMatrixBackend<math::backends::ublas::regular_matrix<int>>);
static_assert(math::ResizableMatrixBackend<math::backends::ublas::compressed_matrix<int>>);
static_assert(math::ResizableVectorBackend<math::backends::ublas::regular_vector<int>>);
static_assert(math::ResizableVectorBackend<math::backends::ublas::compressed_vector<int>>);
static_assert(math::MatrixBackend<math::backends::ublas::regular_matrix<bn254_fp12_value>>);
static_assert(math::MatrixBackend<math::backends::ublas::compressed_matrix<bn254_fp12_value>>);
static_assert(LvalueAddable<math::regular_matrix<int>, math::regular_matrix<int>>);
static_assert(!RvalueAddable<math::regular_matrix<int>, math::regular_matrix<int>>);
static_assert(LvalueAddable<math::regular_vector<int>, math::regular_vector<int>>);
static_assert(!RvalueAddable<math::regular_vector<int>, math::regular_vector<int>>);

BOOST_AUTO_TEST_SUITE(matrix_test_suite)

BOOST_AUTO_TEST_CASE(regular_matrix_access_and_resize) {
    math::regular_matrix<int> value(2, 3);
    BOOST_CHECK_EQUAL(value.rows(), 2);
    BOOST_CHECK_EQUAL(value.columns(), 3);

    value(1, 2) = 7;
    BOOST_CHECK_EQUAL(value(1, 2), 7);

    value.resize(3, 4);
    BOOST_CHECK_EQUAL(value.rows(), 3);
    BOOST_CHECK_EQUAL(value.columns(), 4);
}

BOOST_AUTO_TEST_CASE(compressed_matrix_access_uses_backend_proxy) {
    math::compressed_matrix<int> value(4, 5);
    value(2, 3) = 11;

    BOOST_CHECK_EQUAL(value(2, 3), 11);
    BOOST_CHECK_EQUAL(value(0, 0), 0);
    BOOST_CHECK_EQUAL(value.backend().nnz(), 1);
}

BOOST_AUTO_TEST_CASE(matrix_arithmetic_remains_lazy_until_materialized) {
    math::regular_matrix<int> left(2, 2);
    math::regular_matrix<int> right(2, 2);
    left(0, 0) = 1;
    right(0, 0) = 2;

    auto expression = left + right;
    static_assert(math::MatrixExpression<decltype(expression)>);
    static_assert(!math::MatrixBackend<decltype(expression)>);

    // uBLAS expressions retain references to their operands.
    left(0, 0) = 5;
    math::regular_matrix<int> result(expression);
    BOOST_CHECK_EQUAL(result(0, 0), 7);

    result = left - right;
    BOOST_CHECK_EQUAL(result(0, 0), 3);
}

BOOST_AUTO_TEST_CASE(regular_and_compressed_vectors) {
    math::regular_vector<int> left(3);
    math::regular_vector<int> right(3);
    left[1] = 4;
    right[1] = 6;

    auto expression = left + right;
    static_assert(math::VectorExpression<decltype(expression)>);
    math::regular_vector<int> result(expression);
    BOOST_CHECK_EQUAL(result[1], 10);

    math::compressed_vector<int> sparse(3);
    sparse[2] = 9;
    BOOST_CHECK_EQUAL(sparse[2], 9);
    BOOST_CHECK_EQUAL(sparse.backend().nnz(), 1);
}

BOOST_AUTO_TEST_CASE(bn254_fp12_regular_and_compressed_matrices) {
    const bn254_fp12_value one = bn254_fp12_value::one();
    const bn254_fp12_value two = one.doubled();

    math::regular_matrix<bn254_fp12_value> left(2, 2);
    math::regular_matrix<bn254_fp12_value> right(2, 2);
    left(0, 1) = one;
    right(0, 1) = two;

    auto expression = left + right;
    static_assert(math::MatrixExpression<decltype(expression)>);
    static_assert(!math::MatrixBackend<decltype(expression)>);

    // Confirm that fp12 arithmetic is still lazy through the frontend.
    left(0, 1) = two;
    math::regular_matrix<bn254_fp12_value> result(expression);
    BOOST_CHECK(result(0, 1) == two + two);

    math::compressed_matrix<bn254_fp12_value> sparse(2, 2);
    sparse(1, 0) = result(0, 1);
    BOOST_CHECK(sparse(1, 0) == two + two);
    BOOST_CHECK(sparse(0, 0) == bn254_fp12_value::zero());
    BOOST_CHECK_EQUAL(sparse.backend().nnz(), 1);
}

BOOST_AUTO_TEST_SUITE_END()
