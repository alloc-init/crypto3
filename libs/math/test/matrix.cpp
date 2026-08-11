//---------------------------------------------------------------------------//
// MIT License
//---------------------------------------------------------------------------//

#define BOOST_TEST_MODULE matrix_test

#include <boost/test/unit_test.hpp>

#include <random>
#include <type_traits>
#include <vector>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/math/matrix/backends/ublas/compressed.hpp>
#include <nil/crypto3/math/matrix/backends/ublas/regular.hpp>
#include <nil/crypto3/math/matrix/random.hpp>
#include <nil/crypto3/math/matrix/solve.hpp>

namespace math = nil::crypto3::math;
namespace fields = nil::crypto3::algebra::fields;

using bn254_fp12_field = fields::fp12_2over3over2<fields::alt_bn128<254>>;
using bn254_fp12_value = typename bn254_fp12_field::value_type;

template<typename T>
struct test_matrix_backend {
    using value_type = T;
    using size_type = std::size_t;

    test_matrix_backend() = default;
    test_matrix_backend(size_type rows, size_type columns) : rows_(rows), columns_(columns), values_(rows * columns) {
    }

    size_type size1() const {
        return rows_;
    }
    size_type size2() const {
        return columns_;
    }
    void resize(size_type rows, size_type columns) {
        rows_ = rows;
        columns_ = columns;
        values_.resize(rows * columns);
    }
    value_type &operator()(size_type row, size_type column) {
        return values_[row * columns_ + column];
    }
    const value_type &operator()(size_type row, size_type column) const {
        return values_[row * columns_ + column];
    }

private:
    size_type rows_ = 0;
    size_type columns_ = 0;
    std::vector<value_type> values_;
};

template<typename T>
struct test_vector_backend {
    using value_type = T;
    using size_type = std::size_t;

    test_vector_backend() = default;
    explicit test_vector_backend(size_type size) : values_(size) {
    }

    size_type size() const {
        return values_.size();
    }
    void resize(size_type size) {
        values_.resize(size);
    }
    value_type &operator()(size_type index) {
        return values_[index];
    }
    const value_type &operator()(size_type index) const {
        return values_[index];
    }

private:
    std::vector<value_type> values_;
};

template<typename Left, typename Right>
concept LvalueAddable = requires(Left &left, Right &right) { left + right; };

template<typename Left, typename Right>
concept RvalueAddable = requires { Left() + Right(); };

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

BOOST_AUTO_TEST_CASE(products_inner_product_and_compound_assignment) {
    math::regular_matrix<int> left(2, 2);
    math::regular_matrix<int> right(2, 2);
    math::regular_vector<int> input(2);
    left(0, 0) = 1;
    left(1, 1) = 2;
    right(0, 0) = 3;
    right(1, 1) = 4;
    input[0] = 5;
    input[1] = 6;

    math::regular_matrix<int> matrix_result(math::product(left, right));
    BOOST_CHECK_EQUAL(matrix_result(0, 0), 3);
    BOOST_CHECK_EQUAL(matrix_result(1, 1), 8);

    math::regular_vector<int> vector_result(math::product(left, input));
    BOOST_CHECK_EQUAL(vector_result[0], 5);
    BOOST_CHECK_EQUAL(vector_result[1], 12);
    BOOST_CHECK_EQUAL(math::inner_product(input, vector_result), 97);

    left += right;
    BOOST_CHECK_EQUAL(left(0, 0), 4);
    BOOST_CHECK_EQUAL(left(1, 1), 6);
}

BOOST_AUTO_TEST_CASE(views_row_swap_identity_and_sparse_traversal) {
    auto value = math::identity_matrix<int>(3);
    value(0, 2) = 7;
    math::swap_rows(value, 0, 1);
    BOOST_CHECK_EQUAL(value(1, 2), 7);

    auto block = math::submatrix(value, 0, 2, 1, 3);
    BOOST_CHECK_EQUAL(math::rows(block), 2);
    BOOST_CHECK_EQUAL(math::columns(block), 2);

    auto last_column = math::column(value, 2);
    auto prefix = math::subvector(last_column, 0, 2);
    BOOST_CHECK_EQUAL(prefix.size(), 2);

    math::compressed_matrix<int> sparse(3, 3);
    sparse(0, 2) = 9;
    sparse(2, 1) = 4;
    BOOST_REQUIRE(math::find_element(sparse, 0, 2));
    BOOST_CHECK_EQUAL(*math::find_element(sparse, 0, 2), 9);

    int sum = 0;
    std::size_t count = 0;
    math::for_each_nonzero(sparse, [&](std::size_t, std::size_t, int element) {
        sum += element;
        ++count;
    });
    BOOST_CHECK_EQUAL(sum, 13);
    BOOST_CHECK_EQUAL(count, 2);
}

BOOST_AUTO_TEST_CASE(regular_matrix_solver) {
    const bn254_fp12_value one = bn254_fp12_value::one();
    const bn254_fp12_value two = one.doubled();

    math::regular_matrix<bn254_fp12_value> coefficients(2, 2);
    coefficients(0, 1) = one;
    coefficients(1, 0) = one;
    coefficients(1, 1) = one;

    math::regular_vector<bn254_fp12_value> expected(2);
    expected[0] = one;
    expected[1] = two;
    math::regular_vector<bn254_fp12_value> right_hand_side(2);
    right_hand_side[0] = two;
    right_hand_side[1] = one + two;

    auto result = math::solve(coefficients, right_hand_side);
    BOOST_REQUIRE(result);
    BOOST_CHECK((*result)[0] == expected[0]);
    BOOST_CHECK((*result)[1] == expected[1]);

    math::regular_matrix<bn254_fp12_value> singular(2, 2);
    BOOST_CHECK(!math::solve(singular, right_hand_side));
}

BOOST_AUTO_TEST_CASE(solver_is_backend_generic) {
    using matrix_type = math::matrix<test_matrix_backend<bn254_fp12_value>>;
    using vector_type = math::vector<test_vector_backend<bn254_fp12_value>>;
    const bn254_fp12_value one = bn254_fp12_value::one();
    const bn254_fp12_value two = one.doubled();

    matrix_type coefficients(2, 2);
    coefficients(0, 1) = one;
    coefficients(1, 0) = one;
    coefficients(1, 1) = one;
    vector_type right_hand_side(2);
    right_hand_side[0] = two;
    right_hand_side[1] = one + two;

    auto result = math::solve(coefficients, right_hand_side);
    BOOST_REQUIRE(result);
    BOOST_CHECK((*result)[0] == one);
    BOOST_CHECK((*result)[1] == two);
}

BOOST_AUTO_TEST_CASE(random_generation_is_frontend_and_engine_generic) {
    using matrix_type = math::matrix<test_matrix_backend<unsigned>>;
    using vector_type = math::vector<test_vector_backend<unsigned>>;
    std::mt19937 engine(17);
    auto sample = [](auto &rng) { return rng() % 100; };

    auto matrix = math::random_matrix<matrix_type>(2, 3, engine, sample);
    auto vector = math::random_vector<vector_type>(4, engine, sample);

    BOOST_CHECK_EQUAL(matrix.rows(), 2);
    BOOST_CHECK_EQUAL(matrix.columns(), 3);
    BOOST_CHECK_EQUAL(vector.size(), 4);
    BOOST_CHECK_LT(matrix(1, 2), 100);
    BOOST_CHECK_LT(vector[3], 100);
}

BOOST_AUTO_TEST_CASE(random_generation_samples_field_elements_by_default) {
    std::mt19937_64 engine(23);
    auto vector = math::random_vector<math::regular_vector<bn254_fp12_value>>(3, engine);
    auto matrix = math::random_matrix<math::compressed_matrix<bn254_fp12_value>>(2, 2, engine);

    BOOST_CHECK_EQUAL(vector.size(), 3);
    BOOST_CHECK_EQUAL(matrix.rows(), 2);
    BOOST_CHECK_EQUAL(matrix.columns(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
