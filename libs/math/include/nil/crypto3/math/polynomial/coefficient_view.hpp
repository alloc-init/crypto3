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

#ifndef CRYPTO3_MATH_POLYNOMIAL_COEFFICIENT_VIEW_HPP
#define CRYPTO3_MATH_POLYNOMIAL_COEFFICIENT_VIEW_HPP

#include <concepts>
#include <cstddef>
#include <span>
#include <type_traits>

#include <nil/crypto3/math/polynomial/concepts.hpp>

namespace nil::crypto3::math {

    /**
     * A read-only, non-owning view of a contiguous coefficient range.
     *
     * A view may describe all coefficients of a polynomial or any contiguous
     * subrange. Creating and slicing a view does not copy coefficients. The
     * viewed storage must outlive the view and its subviews; operations that
     * reallocate that storage invalidate them.
     *
     * A coefficient range is not necessarily a canonical polynomial: it may be
     * empty or end in zero coefficients. Consequently, this type deliberately
     * does not provide polynomial degree or representation metadata.
     *
     * first() and last() require count <= size(). subview() requires
     * offset <= size() and either count == std::dynamic_extent or
     * count <= size() - offset.
     */
    template<typename ValueType>
    class coefficient_view {
    public:
        using value_type = ValueType;
        using size_type = std::size_t;
        using const_reference = const value_type &;
        using const_pointer = const value_type *;
        using const_iterator = typename std::span<const value_type>::iterator;

        constexpr coefficient_view() noexcept = default;

        constexpr coefficient_view(const_pointer data, size_type size) noexcept : coefficients_(data, size) {
        }

        template<typename ElementType, std::size_t Extent>
            requires std::same_as<std::remove_const_t<ElementType>, value_type>
        constexpr coefficient_view(std::span<ElementType, Extent> coefficients) noexcept : coefficients_(coefficients) {
        }

        template<CoefficientPolynomial Polynomial>
            requires std::same_as<typename std::remove_cvref_t<Polynomial>::value_type, value_type> &&
                     requires(Polynomial &polynomial) {
                         { polynomial.data() } -> std::convertible_to<const_pointer>;
                     }
        constexpr explicit coefficient_view(Polynomial &polynomial) noexcept :
            coefficients_(polynomial.data(), polynomial.size()) {
        }

        constexpr const_iterator begin() const noexcept {
            return coefficients_.begin();
        }

        constexpr const_iterator end() const noexcept {
            return coefficients_.end();
        }

        constexpr const_reference operator[](size_type index) const noexcept {
            return coefficients_[index];
        }

        constexpr const_pointer data() const noexcept {
            return coefficients_.data();
        }

        constexpr size_type size() const noexcept {
            return coefficients_.size();
        }

        constexpr bool empty() const noexcept {
            return coefficients_.empty();
        }

        constexpr coefficient_view first(size_type count) const {
            return coefficients_.first(count);
        }

        constexpr coefficient_view last(size_type count) const {
            return coefficients_.last(count);
        }

        constexpr coefficient_view subview(size_type offset, size_type count = std::dynamic_extent) const {
            return coefficients_.subspan(offset, count);
        }

    private:
        std::span<const value_type> coefficients_;
    };

    template<CoefficientPolynomial Polynomial>
    coefficient_view(Polynomial &) -> coefficient_view<typename std::remove_cvref_t<Polynomial>::value_type>;

    template<typename ValueType, std::size_t Extent>
    coefficient_view(std::span<ValueType, Extent>) -> coefficient_view<std::remove_const_t<ValueType>>;

    template<typename ValueType>
    coefficient_view(const ValueType *, std::size_t) -> coefficient_view<ValueType>;

}    // namespace nil::crypto3::math

#endif    // CRYPTO3_MATH_POLYNOMIAL_COEFFICIENT_VIEW_HPP
