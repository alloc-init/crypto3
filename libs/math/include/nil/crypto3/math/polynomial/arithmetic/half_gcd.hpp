//---------------------------------------------------------------------------//
// Copyright (c) 2026
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

#ifndef CRYPTO3_MATH_POLYNOMIAL_HALF_GCD_HPP
#define CRYPTO3_MATH_POLYNOMIAL_HALF_GCD_HPP

#include <cstddef>
#include <stdexcept>
#include <utility>

#include <nil/crypto3/math/polynomial/arithmetic/polynomial_division.hpp>
#include <nil/crypto3/math/polynomial/operations/shift.hpp>

namespace nil::crypto3::math::detail {

    /*
     * Let A and B be the input polynomials, with A larger than B. One Euclidean division writes A = q * B + R,
     * where q is the quotient and R is the remainder, and replaces the pair (A, B) with (B, R). Acting on a column
     * pair, that step is the matrix transformation
     *
     *     [ 0   1 ] [ A ]   [ B         ]
     *     [ 1  -q ] [ B ] = [ A - q * B ].
     *
     * The determinant of each step matrix is -1, so the step is invertible and preserves the GCD. Multiplying the
     * step matrices in execution order produces one matrix that represents several consecutive Euclidean divisions.
     *
     * Half-GCD recursively computes such a matrix from the high coefficient halves of A and B, then applies it to the
     * complete inputs to restore the omitted low coefficients exactly. If this first reduction has not made the
     * second polynomial at most half the original size of A, the algorithm performs one complete Euclidean division,
     * recursively reduces the high halves of the new pair, and composes the two recursive matrices with the step
     * matrix. Small inputs construct the same transformation iteratively as the recursion base case.
     */

    /** The four polynomial entries of an accumulated half-GCD transformation matrix. */
    template<SupportsDivrem Backend>
    struct half_gcd_matrix {
        using polynomial_type = typename Backend::polynomial_type;

        polynomial_type entry_00;
        polynomial_type entry_01;
        polynomial_type entry_10;
        polynomial_type entry_11;
    };

    template<SupportsDivrem Backend>
    void set_half_gcd_identity(half_gcd_matrix<Backend> &matrix) {
        using value_type = typename Backend::polynomial_type::value_type;

        matrix.entry_00.resize(1);
        matrix.entry_00[0] = value_type::one();
        matrix.entry_01.resize(1);
        matrix.entry_01[0] = value_type {};
        matrix.entry_10.resize(1);
        matrix.entry_10[0] = value_type {};
        matrix.entry_11.resize(1);
        matrix.entry_11[0] = value_type::one();
    }

    /**
     * Perform one backend-aware Euclidean division without precomputing an inverse for a basecase division.
     *
     * @throws std::invalid_argument if the divisor is zero or has more coefficients than the dividend.
     */
    template<SupportsDivrem Backend>
    void gcd_divrem_step(typename Backend::polynomial_type &quotient,
                         typename Backend::polynomial_type &remainder_output,
                         const typename Backend::polynomial_type &dividend,
                         const typename Backend::polynomial_type &divisor,
                         polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        if (is_zero(divisor)) {
            throw std::invalid_argument("gcd_divrem_step: divisor must be nonzero");
        }
        if (dividend.size() < divisor.size()) {
            throw std::invalid_argument("gcd_divrem_step: dividend must not be smaller than divisor");
        }

        const std::size_t quotient_size = dividend.size() - divisor.size() + 1;
        if (use_basecase_division(arithmetic_context.options(), divisor.size(), quotient_size)) {
            division(quotient, remainder_output, dividend, divisor);
        } else {
            polynomial_divisor_context<Backend> divisor_context(divisor, quotient_size, arithmetic_context);
            divrem(quotient, remainder_output, dividend, divisor_context, arithmetic_context);
        }
    }

    /** Apply a 2-by-2 polynomial matrix to a polynomial pair. */
    template<SupportsDivrem Backend>
    void apply_half_gcd_matrix(typename Backend::polynomial_type &first_output,
                               typename Backend::polynomial_type &second_output,
                               const half_gcd_matrix<Backend> &matrix,
                               const typename Backend::polynomial_type &first,
                               const typename Backend::polynomial_type &second,
                               polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;

        polynomial_type product_00;
        polynomial_type product_01;
        polynomial_type product_10;
        polynomial_type product_11;
        arithmetic_context.multiply(product_00, matrix.entry_00, first);
        arithmetic_context.multiply(product_01, matrix.entry_01, second);
        arithmetic_context.multiply(product_10, matrix.entry_10, first);
        arithmetic_context.multiply(product_11, matrix.entry_11, second);
        addition(first_output, product_00, product_01);
        addition(second_output, product_10, product_11);
    }

    /**
     * Reconstruct the matrix transformation of a complete pair from its low parts and the already transformed high
     * parts. For A = A_low + X^split * A_high and B = B_low + X^split * B_high, this computes
     *
     *     M(A, B) = M(A_low, B_low) + X^split * M(A_high, B_high).
     */
    template<SupportsDivrem Backend>
    void reconstruct_half_gcd_pair(typename Backend::polynomial_type &first_output,
                                   typename Backend::polynomial_type &second_output,
                                   const half_gcd_matrix<Backend> &matrix,
                                   const typename Backend::polynomial_type &first,
                                   const typename Backend::polynomial_type &second,
                                   const typename Backend::polynomial_type &transformed_high_first,
                                   const typename Backend::polynomial_type &transformed_high_second,
                                   std::size_t split,
                                   polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;

        polynomial_type low_first(first);
        polynomial_type low_second(second);
        truncate(low_first, split);
        truncate(low_second, split);

        polynomial_type transformed_low_first;
        polynomial_type transformed_low_second;
        apply_half_gcd_matrix(transformed_low_first, transformed_low_second, matrix, low_first, low_second,
                              arithmetic_context);

        polynomial_type shifted_high_first;
        polynomial_type shifted_high_second;
        shift_left(shifted_high_first, transformed_high_first, split);
        shift_left(shifted_high_second, transformed_high_second, split);
        addition(first_output, transformed_low_first, shifted_high_first);
        addition(second_output, transformed_low_second, shifted_high_second);
    }

    /**
     * Compose two polynomial-pair transformations as output = left * right. When the result is applied to a pair,
     * right acts first and left acts second.
     */
    template<SupportsDivrem Backend>
    void multiply_half_gcd_matrices(half_gcd_matrix<Backend> &output, const half_gcd_matrix<Backend> &left,
                                    const half_gcd_matrix<Backend> &right,
                                    polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;

        half_gcd_matrix<Backend> result;
        polynomial_type first_product;
        polynomial_type second_product;

        arithmetic_context.multiply(first_product, left.entry_00, right.entry_00);
        arithmetic_context.multiply(second_product, left.entry_01, right.entry_10);
        addition(result.entry_00, first_product, second_product);

        arithmetic_context.multiply(first_product, left.entry_00, right.entry_01);
        arithmetic_context.multiply(second_product, left.entry_01, right.entry_11);
        addition(result.entry_01, first_product, second_product);

        arithmetic_context.multiply(first_product, left.entry_10, right.entry_00);
        arithmetic_context.multiply(second_product, left.entry_11, right.entry_10);
        addition(result.entry_10, first_product, second_product);

        arithmetic_context.multiply(first_product, left.entry_10, right.entry_01);
        arithmetic_context.multiply(second_product, left.entry_11, right.entry_11);
        addition(result.entry_11, first_product, second_product);

        output = std::move(result);
    }

    /**
     * Prepend the Euclidean transformation taking the pair (A, B) to (B, A - quotient * B). The first matrix row
     * becomes the old second row; the new second row is the old first row minus quotient times the old second row.
     */
    template<SupportsDivrem Backend>
    void prepend_euclidean_step(half_gcd_matrix<Backend> &matrix,
                                const typename Backend::polynomial_type &quotient,
                                polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;

        polynomial_type quotient_times_entry_10;
        polynomial_type quotient_times_entry_11;
        polynomial_type new_entry_10;
        polynomial_type new_entry_11;
        arithmetic_context.multiply(quotient_times_entry_10, quotient, matrix.entry_10);
        arithmetic_context.multiply(quotient_times_entry_11, quotient, matrix.entry_11);
        subtraction(new_entry_10, matrix.entry_00, quotient_times_entry_10);
        subtraction(new_entry_11, matrix.entry_01, quotient_times_entry_11);

        matrix.entry_00 = std::move(matrix.entry_10);
        matrix.entry_01 = std::move(matrix.entry_11);
        matrix.entry_10 = std::move(new_entry_10);
        matrix.entry_11 = std::move(new_entry_11);
    }

    /**
     * Perform the half-GCD reduction iteratively for a small subproblem. Starting from the input pair and the identity
     * transformation, apply Euclidean quotient steps until the second polynomial is zero or has at most half as many
     * coefficients as the original first polynomial. When requested, matrix accumulates the same quotient steps.
     */
    template<SupportsDivrem Backend>
    void half_gcd_basecase(half_gcd_matrix<Backend> *matrix, typename Backend::polynomial_type &first_output,
                           typename Backend::polynomial_type &second_output,
                           const typename Backend::polynomial_type &first,
                           const typename Backend::polynomial_type &second,
                           polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;

        const std::size_t target_second_size = first.size() / 2;
        if (matrix != nullptr) {
            set_half_gcd_identity(*matrix);
        }
        first_output = first;
        second_output = second;

        polynomial_type quotient;
        polynomial_type next_remainder;
        while (!is_zero(second_output) && second_output.size() > target_second_size) {
            gcd_divrem_step(quotient, next_remainder, first_output, second_output, arithmetic_context);
            first_output = std::move(second_output);
            second_output = std::move(next_remainder);
            if (matrix != nullptr) {
                prepend_euclidean_step(*matrix, quotient, arithmetic_context);
            }
        }
    }

    /**
     * Recursive implementation of the half-GCD algorithm described above. The entry-point overloads validate the
     * input pair; recursive calls preserve those invariants. The transformation matrix is optional because plain GCD
     * needs only the reduced pair, while a future XGCD will also need the accumulated transformation.
     */
    template<SupportsDivrem Backend>
    void half_gcd_reduce_impl(half_gcd_matrix<Backend> *matrix, typename Backend::polynomial_type &first_output,
                              typename Backend::polynomial_type &second_output,
                              const typename Backend::polynomial_type &first,
                              const typename Backend::polynomial_type &second, std::size_t basecase_cutoff,
                              polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        using polynomial_type = typename Backend::polynomial_type;

        const std::size_t split = first.size() / 2;
        if (is_zero(second) || second.size() <= split) {
            if (matrix != nullptr) {
                set_half_gcd_identity(*matrix);
            }
            first_output = first;
            second_output = second;
            return;
        }
        if (basecase_cutoff != 0 && first.size() <= basecase_cutoff) {
            half_gcd_basecase(matrix, first_output, second_output, first, second, arithmetic_context);
            return;
        }

        polynomial_type high_first;
        polynomial_type high_second;
        shift_right(high_first, first, split);
        shift_right(high_second, second, split);

        half_gcd_matrix<Backend> first_matrix;
        polynomial_type transformed_high_first;
        polynomial_type transformed_high_second;
        half_gcd_reduce_impl(&first_matrix, transformed_high_first, transformed_high_second, high_first, high_second,
                             basecase_cutoff, arithmetic_context);

        polynomial_type reduced_first;
        polynomial_type reduced_second;
        reconstruct_half_gcd_pair(reduced_first, reduced_second, first_matrix, first, second, transformed_high_first,
                                  transformed_high_second, split, arithmetic_context);
        if (is_zero(reduced_second) || reduced_second.size() <= split) {
            if (matrix != nullptr) {
                *matrix = std::move(first_matrix);
            }
            first_output = std::move(reduced_first);
            second_output = std::move(reduced_second);
            return;
        }

        polynomial_type quotient;
        polynomial_type next_remainder;
        gcd_divrem_step(quotient, next_remainder, reduced_first, reduced_second, arithmetic_context);
        prepend_euclidean_step(first_matrix, quotient, arithmetic_context);
        if (is_zero(next_remainder)) {
            if (matrix != nullptr) {
                *matrix = std::move(first_matrix);
            }
            first_output = std::move(reduced_second);
            second_output = std::move(next_remainder);
            return;
        }

        // Let l be reduced_second.size(), m be split, and k be second_split. The second recursive input has l - k
        // coefficients, so that reduction leaves at most floor((l - k) / 2) high coefficients. Reconstruction may
        // restore k low coefficients. Choosing k = 2m - l + 1 makes k + floor((l - k) / 2) equal m, the target size.
        const std::size_t second_split = 2 * split - reduced_second.size() + 1;
        polynomial_type high_reduced_second;
        polynomial_type high_next_remainder;
        shift_right(high_reduced_second, reduced_second, second_split);
        shift_right(high_next_remainder, next_remainder, second_split);

        half_gcd_matrix<Backend> second_matrix;
        half_gcd_reduce_impl(&second_matrix, transformed_high_first, transformed_high_second, high_reduced_second,
                             high_next_remainder, basecase_cutoff, arithmetic_context);

        reconstruct_half_gcd_pair(first_output, second_output, second_matrix, reduced_second, next_remainder,
                                  transformed_high_first, transformed_high_second, second_split, arithmetic_context);
        if (matrix != nullptr) {
            multiply_half_gcd_matrices(*matrix, second_matrix, first_matrix, arithmetic_context);
        }
    }

    template<SupportsDivrem Backend>
    void validate_half_gcd_inputs(const typename Backend::polynomial_type &first,
                                  const typename Backend::polynomial_type &second) {
        using value_type = typename Backend::polynomial_type::value_type;

        const auto is_canonical = [](const typename Backend::polynomial_type &polynomial) {
            return polynomial.size() != 0 &&
                   (polynomial.size() == 1 || polynomial[polynomial.size() - 1] != value_type {});
        };
        if (!is_canonical(first) || !is_canonical(second)) {
            throw std::invalid_argument("half_gcd_reduce: inputs must be canonical polynomials");
        }
        if (is_zero(first)) {
            throw std::invalid_argument("half_gcd_reduce: first polynomial must be nonzero");
        }
        if (!is_zero(second) && first.size() <= second.size()) {
            throw std::invalid_argument("half_gcd_reduce: first polynomial must be larger than second polynomial");
        }
    }

    /**
     * Reduce the input pair and return its accumulated transformation matrix. On return,
     *
     *     [ first_output  ]          [ first  ]
     *     [ second_output ] = matrix [ second ].
     *
     * A zero basecase_cutoff disables the iterative base case.
     */
    template<SupportsDivrem Backend>
    void half_gcd_reduce(half_gcd_matrix<Backend> &matrix, typename Backend::polynomial_type &first_output,
                         typename Backend::polynomial_type &second_output,
                         const typename Backend::polynomial_type &first,
                         const typename Backend::polynomial_type &second, std::size_t basecase_cutoff,
                         polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        validate_half_gcd_inputs<Backend>(first, second);
        half_gcd_reduce_impl(&matrix, first_output, second_output, first, second, basecase_cutoff, arithmetic_context);
    }

    /**
     * Reduce the input pair without returning its transformation matrix. Plain GCD uses this overload to avoid the
     * top-level matrix composition. Recursive subproblems still construct the matrices needed for reconstruction.
     * A zero basecase_cutoff disables the iterative base case.
     */
    template<SupportsDivrem Backend>
    void half_gcd_reduce(typename Backend::polynomial_type &first_output,
                         typename Backend::polynomial_type &second_output,
                         const typename Backend::polynomial_type &first,
                         const typename Backend::polynomial_type &second, std::size_t basecase_cutoff,
                         polynomial_arithmetic::polynomial_context<Backend> &arithmetic_context) {
        validate_half_gcd_inputs<Backend>(first, second);
        half_gcd_reduce_impl<Backend>(nullptr, first_output, second_output, first, second, basecase_cutoff,
                                      arithmetic_context);
    }

}    // namespace nil::crypto3::math::detail

#endif    // CRYPTO3_MATH_POLYNOMIAL_HALF_GCD_HPP
