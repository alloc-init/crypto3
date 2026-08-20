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

#ifndef CRYPTO3_MATH_POLYNOMIAL_FACTORIZATION_HPP
#define CRYPTO3_MATH_POLYNOMIAL_FACTORIZATION_HPP

#include <cstddef>
#include <vector>

#include <nil/crypto3/math/polynomial/concepts.hpp>

namespace nil::crypto3::math {

    /** The action requested by a callback after a factorization stage produces a factor. */
    enum class factorization_control { continue_factorization, stop_factorization };

    /** A polynomial factor together with its positive multiplicity in the input polynomial. */
    template<CoefficientPolynomial Polynomial>
    struct polynomial_factor {
        Polynomial polynomial;
        std::size_t multiplicity = 1;

        bool operator==(const polynomial_factor &) const = default;
    };

    /**
     * A collected polynomial factorization. A complete result satisfies
     *
     *     input = leading_coefficient * product(factor.polynomial ^ factor.multiplicity).
     *
     * Nonconstant factors are stored in monic canonical form. When a staged callback requests an early stop, factors
     * contains the produced prefix, including the factor that caused the stop, and complete is false.
     */
    template<CoefficientPolynomial Polynomial>
    struct polynomial_factorization_result {
        using polynomial_type = Polynomial;
        using value_type = typename polynomial_type::value_type;
        using factor_type = polynomial_factor<polynomial_type>;

        value_type leading_coefficient {};
        std::vector<factor_type> factors;
        bool complete = true;

        bool operator==(const polynomial_factorization_result &) const = default;
    };

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_FACTORIZATION_HPP
