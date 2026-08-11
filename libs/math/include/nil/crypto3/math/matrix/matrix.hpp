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
    template<MatrixBackend Backend>
    class matrix {
    public:
        using backend_type = Backend;
        using value_type = typename backend_type::value_type;
        using size_type = typename backend_type::size_type;

        matrix() = default;

        matrix(size_type rows, size_type columns)
            requires std::constructible_from<backend_type, size_type, size_type>
            : backend_(rows, columns) {
        }

        explicit matrix(backend_type backend) : backend_(std::move(backend)) {
        }

        // This is the materialization boundary for a lazy uBLAS expression.
        template<MatrixExpression Expression>
            requires(!std::same_as<std::remove_cvref_t<Expression>, matrix> &&
                     std::constructible_from<backend_type, Expression &&>)
        matrix(Expression &&expression) : backend_(std::forward<Expression>(expression)) {
        }

        template<MatrixExpression Expression>
            requires requires(backend_type &backend, Expression &&expression) {
                backend = std::forward<Expression>(expression);
            }
        matrix &operator=(Expression &&expression) {
            backend_ = std::forward<Expression>(expression);
            return *this;
        }

        size_type rows() const noexcept {
            return backend_.size1();
        }

        size_type columns() const noexcept {
            return backend_.size2();
        }

        void resize(size_type rows, size_type columns)
            requires ResizableMatrixBackend<backend_type> {
            backend_.resize(rows, columns);
        }

        // decltype(auto) preserves dense references and sparse element proxies.
        decltype(auto) operator()(size_type row, size_type column) {
            return backend_(row, column);
        }

        decltype(auto) operator()(size_type row, size_type column) const {
            return backend_(row, column);
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

    template<MatrixBackend Backend>
    matrix(Backend) -> matrix<Backend>;

}    // namespace nil::crypto3::math
