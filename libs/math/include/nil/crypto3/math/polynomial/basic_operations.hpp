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

#ifndef CRYPTO3_MATH_POLYNOMIAL_BASIC_OPERATIONS_HPP
#define CRYPTO3_MATH_POLYNOMIAL_BASIC_OPERATIONS_HPP

#include <algorithm>
#include <concepts>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/math/algorithms/unity_root.hpp>
#include <nil/crypto3/math/domains/detail/basic_radix2_domain_aux.hpp>
#include <nil/crypto3/math/detail/field_utils.hpp>
#include <nil/crypto3/math/polynomial/polynomial_backend.hpp>

namespace nil {
    namespace crypto3 {
        namespace math {
            namespace detail {
                template<typename Range>
                concept PolynomialCoefficientRange =
                    std::ranges::random_access_range<const Range> && std::ranges::sized_range<const Range>;

                template<typename Range>
                concept MutablePolynomialCoefficientRange =
                    PolynomialCoefficientRange<Range> && std::ranges::random_access_range<Range> &&
                    requires(Range &range, const Range &source, std::ranges::range_size_t<const Range> size,
                             const std::ranges::range_value_t<Range> &value) {
                        range = source;
                        range[size] = value;
                        range.resize(size);
                        range.resize(size, value);
                    };

            }    // namespace detail

            /**
             * Returns true if polynomial A is a zero polynomial.
             */
            template<std::input_iterator Iter>
                requires std::default_initializable<std::iter_value_t<Iter>> &&
                         std::equality_comparable<std::iter_value_t<Iter>>
            bool is_zero(const Iter &begin, const Iter &end) {
                std::iter_value_t<Iter> zero {};
                return !std::any_of(begin, end, [&zero](const auto &value) { return value != zero; });
            }

            /**
             * Returns true if polynomial A is a zero polynomial.
             */
            template<typename Range>
                requires std::ranges::input_range<const Range> &&
                         std::default_initializable<std::ranges::range_value_t<const Range>> &&
                         std::equality_comparable<std::ranges::range_value_t<const Range>>
            bool is_zero(const Range &a) {
                return is_zero(std::ranges::begin(a), std::ranges::end(a));
            }

            /**
             * Reverse the coefficient order, then truncate the result to n coefficients.
             *
             * @pre n <= a.size().
             */
            template<detail::MutablePolynomialCoefficientRange Range>
            void reverse(Range &a, std::size_t n) {
                // Legacy std::reverse only requires swappable elements. std::ranges::reverse additionally requires
                // std::permutable, which some supported Crypto3 curve elements do not satisfy.
                std::reverse(std::ranges::begin(a), std::ranges::end(a));
                a.resize(n);
            }

            namespace detail {
                template<MutablePolynomialCoefficientRange AlgebraicRange, typename MultiplyOperation>
                    requires std::copy_constructible<AlgebraicRange> && std::default_initializable<AlgebraicRange> &&
                             requires(AlgebraicRange &result, const AlgebraicRange &input,
                                      MultiplyOperation &multiply_operation) {
                                 multiply_operation(result, input);
                                 {
                                     std::ranges::range_value_t<AlgebraicRange>::zero()
                                 } -> std::convertible_to<std::ranges::range_value_t<AlgebraicRange>>;
                             }
                AlgebraicRange transpose_multiplication_impl(const std::size_t n, const AlgebraicRange &a,
                                                             MultiplyOperation multiply_operation) {
                    const std::size_t m = std::ranges::size(a);

                    AlgebraicRange product(a);
                    reverse(product, m);
                    multiply_operation(product, product);
                    product.resize(m + n, std::ranges::range_value_t<AlgebraicRange>::zero());

                    AlgebraicRange result;
                    result.resize(n + 1, std::ranges::range_value_t<AlgebraicRange>::zero());
                    for (std::size_t i = 0; i <= n; ++i) {
                        result[i] = product[m - 1 + i];
                    }
                    return result;
                }
            }    // namespace detail

            /**
             * Removes trailing zero coefficients while retaining at least one coefficient.
             *
             * Example - Degree-4 Polynomial: [0, 1, 2, 3, 4, 0, 0, 0, 0] -> [0, 1, 2, 3, 4]
             * Empty inputs and all representations of the zero polynomial become [0].
             */
            template<detail::MutablePolynomialCoefficientRange Range>
                requires std::default_initializable<std::ranges::range_value_t<Range>> &&
                         std::equality_comparable<std::ranges::range_value_t<Range>>
            void condense(Range &a) {
                std::size_t i = std::ranges::size(a);
                std::ranges::range_value_t<Range> zero {};
                if (i == 0) {
                    a.resize(1, zero);
                    return;
                }
                while (i > 1 && a[i - 1] == zero) {
                    --i;
                }
                a.resize(i);
            }

            /**
             * Replace a coefficient polynomial by its first coefficient_count coefficients and normalize the
             * result. This computes the polynomial modulo X^coefficient_count; a zero coefficient_count produces
             * [0].
             */
            template<CoefficientPolynomial Polynomial>
                requires detail::MutablePolynomialCoefficientRange<Polynomial> &&
                         std::default_initializable<typename Polynomial::value_type> &&
                         std::equality_comparable<typename Polynomial::value_type>
            void truncate(Polynomial &polynomial, std::size_t coefficient_count) {
                using value_type = typename Polynomial::value_type;

                if (coefficient_count == 0) {
                    polynomial.resize(1);
                    polynomial[0] = value_type {};
                    return;
                }
                if (polynomial.size() > coefficient_count) {
                    polynomial.resize(coefficient_count);
                }
                condense(polynomial);
            }

            /**
             * Computes the standard polynomial addition, polynomial A + polynomial B, and stores result in
             * polynomial C.
             * The output may alias either input.
             *
             * @pre a and b are not empty.
             */
            template<detail::MutablePolynomialCoefficientRange Range>
                requires std::default_initializable<std::ranges::range_value_t<Range>> &&
                         std::equality_comparable<std::ranges::range_value_t<Range>> &&
                         requires(const std::ranges::range_value_t<Range> &left,
                                  const std::ranges::range_value_t<Range> &right) {
                             { left + right } -> std::convertible_to<std::ranges::range_value_t<Range>>;
                         }
            void addition(Range &c, const Range &a, const Range &b) {

                using value_type = std::ranges::range_value_t<Range>;

                if (is_zero(a)) {
                    c = b;
                } else if (is_zero(b)) {
                    c = a;
                } else {
                    std::size_t a_size = std::ranges::size(a);
                    std::size_t b_size = std::ranges::size(b);

                    if (a_size > b_size) {
                        c.resize(a_size);
                        std::transform(std::ranges::begin(b), std::ranges::end(b), std::ranges::begin(a),
                                       std::ranges::begin(c), std::plus<value_type>());
                        std::copy(std::ranges::begin(a) + b_size, std::ranges::end(a), std::ranges::begin(c) + b_size);
                    } else {
                        c.resize(b_size);
                        std::transform(std::ranges::begin(a), std::ranges::end(a), std::ranges::begin(b),
                                       std::ranges::begin(c), std::plus<value_type>());
                        std::copy(std::ranges::begin(b) + a_size, std::ranges::end(b), std::ranges::begin(c) + a_size);
                    }
                }

                condense(c);
            }

            /**
             * Computes the standard polynomial subtraction, polynomial A - polynomial B, and stores result in
             * polynomial C.
             * The output may alias either input.
             *
             * @pre a and b are not empty.
             */
            template<detail::MutablePolynomialCoefficientRange Range>
                requires std::default_initializable<std::ranges::range_value_t<Range>> &&
                         std::equality_comparable<std::ranges::range_value_t<Range>> &&
                         requires(const std::ranges::range_value_t<Range> &left,
                                  const std::ranges::range_value_t<Range> &right) {
                             { left - right } -> std::convertible_to<std::ranges::range_value_t<Range>>;
                             { -right } -> std::convertible_to<std::ranges::range_value_t<Range>>;
                         }
            void subtraction(Range &c, const Range &a, const Range &b) {

                using value_type = std::ranges::range_value_t<Range>;

                if (is_zero(b)) {
                    c = a;
                } else if (is_zero(a)) {
                    c.resize(std::ranges::size(b));
                    std::transform(std::ranges::begin(b), std::ranges::end(b), std::ranges::begin(c),
                                   std::negate<value_type>());
                } else {
                    std::size_t a_size = std::ranges::size(a);
                    std::size_t b_size = std::ranges::size(b);

                    if (a_size > b_size) {
                        c.resize(a_size);
                        std::transform(std::ranges::begin(a), std::ranges::begin(a) + b_size, std::ranges::begin(b),
                                       std::ranges::begin(c), std::minus<value_type>());
                        std::copy(std::ranges::begin(a) + b_size, std::ranges::end(a), std::ranges::begin(c) + b_size);
                    } else {
                        c.resize(b_size);
                        std::transform(std::ranges::begin(a), std::ranges::end(a), std::ranges::begin(b),
                                       std::ranges::begin(c), std::minus<value_type>());
                        std::transform(std::ranges::begin(b) + a_size, std::ranges::end(b),
                                       std::ranges::begin(c) + a_size, std::negate<value_type>());
                    }
                }

                condense(c);
            }

            /**
             * Multiply two polynomials using a radix-2 FFT whose size is the smallest power of two containing the
             * complete product, and store the result in polynomial C.
             * This handles both ordinary field-by-field polynomial multiplication and curve-by-field module-valued
             * convolution. FieldRange contains field elements, while AlgebraicRange contains either field or curve
             * elements.
             * The output may alias either input when their range types permit it.
             *
             * @pre a and b are not empty.
             * @todo Move ordinary same-coefficient-type multiplication into a radix-2 backend. Preserve the
             * curve-by-field case separately as module-valued convolution.
             */
            template<detail::MutablePolynomialCoefficientRange AlgebraicRange,
                     detail::MutablePolynomialCoefficientRange FieldRange>
                requires requires { typename std::ranges::range_value_t<FieldRange>::field_type; } &&
                         algebra::is_field<typename std::ranges::range_value_t<FieldRange>::field_type>::value &&
                         std::same_as<typename std::ranges::range_value_t<FieldRange>::field_type::value_type,
                                      std::ranges::range_value_t<FieldRange>> &&
                         std::default_initializable<std::ranges::range_value_t<AlgebraicRange>> &&
                         std::equality_comparable<std::ranges::range_value_t<AlgebraicRange>> &&
                         requires(std::ranges::range_value_t<AlgebraicRange> &algebraic_value,
                                  const std::ranges::range_value_t<AlgebraicRange> &const_algebraic_value,
                                  const std::ranges::range_value_t<FieldRange> &field_value) {
                             {
                                 std::ranges::range_value_t<AlgebraicRange>::zero()
                             } -> std::convertible_to<std::ranges::range_value_t<AlgebraicRange>>;
                             {
                                 const_algebraic_value * field_value
                             } -> std::convertible_to<std::ranges::range_value_t<AlgebraicRange>>;
                             algebraic_value *= field_value;
                             algebraic_value += const_algebraic_value;
                             algebraic_value -= const_algebraic_value;
                         }
            void multiplication(AlgebraicRange &c, const AlgebraicRange &a, const FieldRange &b) {
                using algebraic_value_type = std::ranges::range_value_t<AlgebraicRange>;
                using field_value_type = std::ranges::range_value_t<FieldRange>;

                typedef typename field_value_type::field_type FieldType;
                BOOST_ASSERT_MSG(!std::ranges::empty(a), "Uninitialized polynomial");
                BOOST_ASSERT_MSG(!std::ranges::empty(b), "Uninitialized polynomial");

                const std::size_t n = detail::power_of_two(std::ranges::size(a) + std::ranges::size(b) - 1);
                field_value_type omega = unity_root<FieldType>(n);

                AlgebraicRange u(a);
                FieldRange v(b);
                u.resize(n, algebraic_value_type::zero());
                v.resize(n, field_value_type::zero());
                c.resize(n, algebraic_value_type::zero());

                detail::basic_radix2_fft<FieldType>(u, omega);
                detail::basic_radix2_fft<FieldType>(v, omega);

                for (std::size_t i = 0; i < n; ++i) {
                    c[i] = u[i] * v[i];
                }

                detail::basic_radix2_fft<FieldType>(c, omega.inversed());

                const field_value_type size_inverse = field_value_type(n).inversed();

                for (std::size_t i = 0; i < n; ++i) {
                    c[i] = c[i] * size_inverse;
                }

                condense(c);
            }

            /**
             * Multiply canonical coefficient polynomials using a reusable polynomial-arithmetic context.
             * The output is canonical and may alias either input.
             */
            template<polynomial_arithmetic::PolynomialBackend Backend>
            void multiplication(typename Backend::polynomial_type &output,
                                const typename Backend::polynomial_type &left,
                                const typename Backend::polynomial_type &right,
                                polynomial_arithmetic::polynomial_context<Backend> &context) {
                context.multiply(output, left, right);
            }

            /**
             * Square a canonical coefficient polynomial using a reusable polynomial-arithmetic context.
             * The output is canonical and may alias the input.
             */
            template<polynomial_arithmetic::PolynomialBackend Backend>
            void square(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &input,
                        polynomial_arithmetic::polynomial_context<Backend> &context) {
                context.square(output, input);
            }

            /**
             * Compute the product modulo X^coefficient_count using a reusable polynomial-arithmetic context.
             * The output is canonical and may alias either input. A coefficient count of zero produces [0].
             */
            template<polynomial_arithmetic::PolynomialBackend Backend>
            void multiply_low(typename Backend::polynomial_type &output, const typename Backend::polynomial_type &left,
                              const typename Backend::polynomial_type &right, std::size_t coefficient_count,
                              polynomial_arithmetic::polynomial_context<Backend> &context) {
                context.multiply_low(output, left, right, coefficient_count);
            }

            /**
             * Compute the transposed multiplication of a by c. If m is a.size(), reverse a, multiply it by c, and
             * return the n + 1 coefficients at product indices m - 1 through m + n - 1.
             *
             * Below we use the transposed multiplication definition from
             * [Bostan, Lecerf, & Schost, 2003. Tellegen's Principle in Practice, on page 39].
             *
             * @pre a is not empty.
             * @note The internal convolution currently uses the legacy radix-2 multiplication overload.
             */
            template<detail::MutablePolynomialCoefficientRange AlgebraicRange,
                     detail::MutablePolynomialCoefficientRange FieldRange>
                requires std::copy_constructible<AlgebraicRange> && std::default_initializable<AlgebraicRange> &&
                         requires(AlgebraicRange &result, const AlgebraicRange &input, const FieldRange &coefficients) {
                             multiplication(result, input, coefficients);
                         }
            AlgebraicRange transpose_multiplication(const std::size_t &n, const AlgebraicRange &a,
                                                    const FieldRange &c) {
                return detail::transpose_multiplication_impl(
                    n, a,
                    [&c](AlgebraicRange &output, const AlgebraicRange &input) { multiplication(output, input, c); });
            }

            /**
             * Compute transposed multiplication using a reusable polynomial-arithmetic context. Field coefficients
             * are embedded into the backend's coefficient ring before multiplication.
             *
             * @pre a is not empty.
             */
            template<polynomial_arithmetic::PolynomialBackend Backend, detail::PolynomialCoefficientRange FieldRange>
                requires std::default_initializable<typename Backend::polynomial_type> &&
                         algebra::is_field_element<std::ranges::range_value_t<const FieldRange>>::value &&
                         requires(const std::ranges::range_value_t<const FieldRange> &field_value) {
                             {
                                 Backend::polynomial_type::value_type::one() * field_value
                             } -> std::convertible_to<typename Backend::polynomial_type::value_type>;
                         }
            typename Backend::polynomial_type
                transpose_multiplication(const std::size_t n, const typename Backend::polynomial_type &a,
                                         const FieldRange &c,
                                         polynomial_arithmetic::polynomial_context<Backend> &context) {
                using polynomial_type = typename Backend::polynomial_type;
                using value_type = typename polynomial_type::value_type;

                polynomial_type embedded_coefficients(std::ranges::size(c), value_type::zero());
                std::size_t coefficient_index = 0;
                for (const auto &coefficient : c) {
                    embedded_coefficients[coefficient_index++] = value_type::one() * coefficient;
                }

                return detail::transpose_multiplication_impl(
                    n, a, [&embedded_coefficients, &context](polynomial_type &output, const polynomial_type &input) {
                        multiplication(output, input, embedded_coefficients, context);
                    });
            }

            /**
             * Perform quadratic Euclidean polynomial division, producing Q and R such that A = Q * B + R.
             * This implementation is suitable as the small-degree fallback for future fast division algorithms.
             * Existing contents of Q and R are replaced.
             *
             * @pre A and B are nonempty canonical coefficient vectors.
             * @pre B has a nonzero leading coefficient.
             * @pre Q and R are distinct and do not alias either input.
             */
            template<detail::MutablePolynomialCoefficientRange Range>
                requires algebra::is_field_element<std::ranges::range_value_t<Range>>::value &&
                         std::copy_constructible<Range> &&
                         std::constructible_from<Range, std::ranges::range_size_t<const Range>,
                                                 std::ranges::range_value_t<Range>> &&
                         std::constructible_from<Range, std::ranges::iterator_t<const Range>,
                                                 std::ranges::sentinel_t<const Range>>
            void division(Range &q, Range &r, const Range &a, const Range &b) {
                using value_type = std::ranges::range_value_t<Range>;

                const std::size_t divisor_degree = std::ranges::size(b) - 1;

                // Special case when B has degree 0.
                if (divisor_degree == 0) {
                    const value_type divisor_inverse = b[0].inversed();
                    q.resize(std::ranges::size(a));
                    std::transform(std::ranges::begin(a), std::ranges::end(a), std::ranges::begin(q),
                                   [&divisor_inverse](const value_type &value) { return value * divisor_inverse; });
                    // Division by a nonzero constant always has zero remainder.
                    r.resize(1);
                    r[0] = value_type::zero();
                }
                // Divide by B(X) = X^N + C in linear time. Modulo B, X^N = -C, so each
                // coefficient above degree N - 1 is moved into the quotient and folded N
                // positions down into the remaining dividend by multiplication with -C.
                else if (b[std::ranges::size(b) - 1] == value_type::one() &&
                         is_zero(std::ranges::begin(b) + 1, std::ranges::end(b) - 1) &&
                         std::ranges::size(a) >= std::ranges::size(b)) {
                    q = Range(std::ranges::size(a) - std::ranges::size(b) + 1, value_type::zero());
                    r = Range(std::ranges::begin(a),
                              std::ranges::end(a) - (std::ranges::size(a) - std::ranges::size(b) + 1));

                    const value_type negated_constant = -b[0];
                    auto end = --std::ranges::end(a);
                    for (std::size_t t = std::ranges::size(q); t != 0; --t, --end) {
                        q[t - 1] += *end;
                        if (t - 1 >= divisor_degree) {
                            q[t - 1 - divisor_degree] = q[t - 1] * negated_constant;
                        } else {
                            r[t - 1] += q[t - 1] * negated_constant;
                        }
                    }
                    condense(r);
                } else {
                    const value_type inverse_leading_coefficient = b[std::ranges::size(b) - 1].inversed();
                    r = Range(a);
                    q = Range(std::ranges::size(r), value_type::zero());

                    std::size_t remainder_degree = std::ranges::size(r) - 1;

                    while (remainder_degree >= divisor_degree && !is_zero(r)) {
                        const std::size_t degree_shift = remainder_degree - divisor_degree;
                        const value_type quotient_coefficient = r[remainder_degree] * inverse_leading_coefficient;

                        q[degree_shift] = quotient_coefficient;

                        std::transform(std::ranges::begin(b), std::ranges::end(b), std::ranges::begin(r) + degree_shift,
                                       std::ranges::begin(r) + degree_shift,
                                       [&quotient_coefficient](const value_type &divisor_coefficient,
                                                               const value_type &remainder_coefficient) {
                                           return remainder_coefficient - divisor_coefficient * quotient_coefficient;
                                       });

                        condense(r);
                        remainder_degree = std::ranges::size(r) - 1;
                    }
                }
                condense(q);
            }
        }    // namespace math
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_MATH_POLYNOMIAL_BASIC_OPERATIONS_HPP
