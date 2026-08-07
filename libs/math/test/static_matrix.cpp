//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//---------------------------------------------------------------------------//

#define BOOST_TEST_MODULE constexpr_matrix_test

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <nil/crypto3/math/static_matrix/matrix.hpp>
#include <nil/crypto3/math/static_matrix/math.hpp>
#include <nil/crypto3/math/static_matrix/operators.hpp>

#include <nil/crypto3/algebra/fields/goldilocks.hpp>

using namespace nil::crypto3::algebra;
using namespace nil::crypto3::math;

using field = fields::goldilocks;
using value = field::value_type;


static_assert(make_vector(1, 2, 3) == vector {1, 2, 3}, "make_vector and uniform initialization deduction guide");

static_assert(make_vector(1, 2, 3) == vector {{1, 2, 3}}, "make_vector and aggregate initialization deduction guide");

static_assert(elementwise([](value x) { return x.inversed(); }, vector {1, 2, 4}) ==
                  vector<value, 3> {1, 0x7fffffff80000001_cppui_modular64, 0xbfffffff40000001_cppui_modular64},
              "elementwise");

static_assert(vector {1, 2, 3} == vector {1, 2, 3}, "operator==");

static_assert(vector {1, 2, 3} != vector {3, 2, 1}, "operator!=");

static_assert(vector {1, 2, 3} + vector {1, 2, 3} == vector {2, 4, 6}, "operator+");

static_assert(sum(vector {1, 2, 3}) == 6, "sum");

static_assert(iota<5>(0) == vector {0, 1, 2, 3, 4}, "iota");

static_assert(iota<5, value>() == vector<value, 5> {0, 1, 2, 3, 4}, "iota");

static_assert(fill<4>(2) == vector {2, 2, 2, 2}, "fill");

static_assert(generate<4>([](auto i) { return value(i * i); }) == vector<value, 4> {0, 1, 4, 9}, "generate");

static_assert(vector {1, 2, 3} == slice<3>(vector {1, 2, 3, 4}), "slice-no offset");

static_assert(vector {2, 3, 4} == slice<3>(vector {1, 2, 3, 4}, 1), "slice with offset");


// Uniform initialization
constexpr matrix<value, 3, 3> m1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};

// Type deduction
constexpr matrix m2 = {{{value(1), value(2)}}};

constexpr matrix m22 = {{{value(1), value(3)}, {value(2), value(7)}}};

static_assert(m1[0][2] == 3, "matrix[]");

static_assert(m1.row(2) == vector {value(7), value(8), value(9)}, "matrix row");

static_assert(m1.column(2) == vector {value(3), value(6), value(9)}, "matrix column");

static_assert(fill<2, 2>(value(3)) == matrix {{{value(3), value(3)}, {value(3), value(3)}}}, "matrix fill");

static_assert(matmul(m1, m1) == matrix<value, 3, 3> {{{30, 36, 42}, {66, 81, 96}, {102, 126, 150}}},
              "real matrix multiply");

static_assert(identity<value, 3> == matrix<value, 3, 3> {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}, "identity");

static_assert(identity<value, 3> == inverse(identity<value, 3>), "inverse-identity");

static_assert(inverse(m22) == matrix<value, 2, 2> {{{7, -3}, {-2, 1}}}, "inverse");

static_assert(matmul(inverse(m22), matrix<value, 2, 1> {{{1}, {1}}}) == matrix<value, 2, 1> {{{4}, {-1}}},
              "A^-1*b = x");

static_assert(horzcat(identity<value, 2>, identity<value, 2>) == matrix<value, 2, 4> {{{1, 0, 1, 0}, {0, 1, 0, 1}}},
              "horzcat");

static_assert(submat<2, 2>(m1, 1, 1) == matrix<value, 2, 2> {{{5, 6}, {8, 9}}}, "submat");

static_assert(rref(m1) == matrix<value, 3, 3> {{{1, 0, -1}, {0, 1, 2}, {0, 0, 0}}}, "rref");

static_assert(rank(m1) == 2, "rank");

