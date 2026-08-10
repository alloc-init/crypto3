//---------------------------------------------------------------------------//
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
    template<typename T, std::size_t N, VectorBackend Backend>
    struct vector {
        static_assert(N != 0, "vector must contain at least one element");

        T &operator[](size_t i) noexcept {
            return backend[i];
        }

        const T &operator[](size_t i) const noexcept {
            return backend[i];
        }

        constexpr auto begin() noexcept {
            return backend.begin();
        }

        constexpr auto end() noexcept {
            return backend.end();
        }

        constexpr auto begin() const noexcept {
            return backend.begin();
        }

        constexpr auto end() const noexcept {
            return backend.end();
        }

        constexpr auto cbegin() const noexcept {
            return backend.cbegin();
        }

        constexpr auto cend() const noexcept {
            return backend.cend();
        }

        Backend backend;
    };

    template<typename... Args>
    constexpr decltype(auto) make_vector(Args... args) {
        return vector {args...};
    }

    template<typename T, typename... U>
        requires(std::same_as<T, U> && ...)
    vector(T, U...) -> vector<T, 1 + sizeof...(U)>;

    template<typename T, std::size_t N>
    vector(const T (&)[N]) -> vector<T, N>;

}    // namespace nil::crypto3::math
