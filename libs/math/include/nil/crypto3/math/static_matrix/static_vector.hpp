//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2020-2021 Ilias Khairullin <ilias@nil.foundation>
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

#pragma once

#include <array>
#include <cstddef>

namespace nil::crypto3::math {
    /** @brief A container representing a vector
     *    @tparam T scalar type to contain
     *    @tparam N size of the vector
     *
     *    `static_vector` is a container representing a vector at compile time.
     *
     *    It is an aggregate type similar to `std::array`, and can be initialized with
     *    aggregate initialization or with the `make_vector` function.
     */
    template<typename T, std::size_t N>
    struct static_vector {
        static_assert(N != 0, "static vector must contain at least one element");

        using value_type = T;
        using size_type = std::size_t;
        static constexpr size_type size = N;    ///< @brief size of the vector

        /** @name Element access */
        ///@{
        /** @brief access specified element
         *    @param i position of the scalar element
         *    @return the requested scalar element
         *
         *    Returns a reference to the scalar element in position `i`, without bounds checking.
         */
        constexpr T &operator[](size_type i) noexcept {
            return array[i];
        }

        /// @copydoc operator[]
        constexpr const T &operator[](size_type i) const noexcept {
            return array[i];
        }

        ///@}

        /** @name Iterators */
        ///@{
        /** @brief returns an iterator to the beginning
         *    @return an iterator to the beginning
         *
         *    Returns an iterator to the beginning of the vector.
         */
        constexpr T *begin() noexcept {
            return array;
        }

        /** @brief returns an iterator to the end
         *    @return an iterator past the end
         *
         *    Returns an iterator to the end of the vector.
         */
        constexpr T *end() noexcept {
            return array + N;
        }

        /// @copydoc begin
        constexpr const T *cbegin() const noexcept {
            return array;
        }

        /// @copydoc end
        constexpr const T *cend() const noexcept {
            return array + N;
        }

        ///@}

        T array[N];    ///< @private
    };

    /** \addtogroup static_vector
     *    @{
     */

    /** @brief constructs a `static_vector` from arguments
     *    @param args scalar elements to combine into a static_vector
     *    @return a static_vector containing `args`
     *    @relatesalso static_vector
     *
     *    Constructs a static_vector from its arguments, checking that all arguments are of
     *    the same type.
     */
    template<typename... Args>
    constexpr decltype(auto) make_vector(Args... args) {
        return static_vector {args...};
    }

    /** @name static_vector deduction guides */
    ///@{

    /** @brief deduction guide for uniform initialization
     *    @relatesalso static_vector
     *
     *    This deduction guide allows static_vector to be constructed like this:
     *    \code{.cpp}
     *    static_vector v{1., 2.}; // deduces the type of v to be static_vector<double, 2>
     *    \endcode
     */
    template<typename T, typename... U>
    static_vector(T, U...)
        -> static_vector<typename std::enable_if<(std::is_same<T, U>::value && ...), T>::type, 1 + sizeof...(U)>;

    /** @brief deduction guide for aggregate initialization
     *    @relatesalso static_vector
     *
     *    This deduction guide allows static_vector to be constructed like this:
     *    \code{.cpp}
     *    static_vector v{{1., 2.}}; // deduces the type of v to be static_vector<double, 2>
     *    \endcode
     */
    template<typename T, std::size_t N>
    static_vector(const T (&)[N]) -> static_vector<T, N>;

    ///@}

    /** @}*/
}    // namespace nil::crypto3::math
