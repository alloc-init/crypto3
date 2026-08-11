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

#include <cstddef>
#include <type_traits>
#include <utility>

#include <nil/crypto3/math/matrix/concepts.hpp>

namespace nil::crypto3::math {
    template<VectorBackend Backend>
    class vector {
    public:
        using backend_type = Backend;
        using value_type = typename backend_type::value_type;
        using size_type = typename backend_type::size_type;

        vector() = default;

        explicit vector(size_type size)
            requires std::constructible_from<backend_type, size_type>
            : backend_(size) {
        }

        explicit vector(backend_type backend) : backend_(std::move(backend)) {
        }

        template<VectorExpression Expression>
            requires(!std::same_as<std::remove_cvref_t<Expression>, vector> &&
                     std::constructible_from<backend_type, Expression &&>)
        vector(Expression &&expression) : backend_(std::forward<Expression>(expression)) {
        }

        template<VectorExpression Expression>
            requires requires(backend_type &backend, Expression &&expression) {
                backend = std::forward<Expression>(expression);
            }
        vector &operator=(Expression &&expression) {
            backend_ = std::forward<Expression>(expression);
            return *this;
        }

        size_type size() const noexcept {
            return backend_.size();
        }

        void resize(size_type size)
            requires ResizableVectorBackend<backend_type> {
            backend_.resize(size);
        }

        decltype(auto) operator[](size_type i) {
            return backend_(i);
        }

        decltype(auto) operator[](size_type i) const {
            return backend_(i);
        }

        backend_type &backend() & noexcept {
            return backend_;
        }

        const backend_type &backend() const & noexcept {
            return backend_;
        }

        backend_type &&backend() && noexcept {
            return std::move(backend_);
        }

    private:
        backend_type backend_;
    };

    template<VectorBackend Backend>
    vector(Backend) -> vector<Backend>;

}    // namespace nil::crypto3::math
