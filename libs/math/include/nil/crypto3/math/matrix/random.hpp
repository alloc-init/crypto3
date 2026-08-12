#pragma once

#include <concepts>
#include <cstddef>
#include <functional>
#include <utility>

#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/math/matrix/matrix.hpp>
#include <nil/crypto3/math/matrix/vector.hpp>

namespace nil::crypto3::math {
    template<typename Vector, typename RandomEngine, typename Sampler>
        requires std::constructible_from<Vector, std::size_t> &&
                 requires(Vector &value, RandomEngine &engine, Sampler &sampler) {
                     value[0] = std::invoke(sampler, engine);
                 }
    Vector random_vector(std::size_t size, RandomEngine &engine, Sampler sampler) {
        Vector result(size);
        for (std::size_t i = 0; i < size; ++i) {
            result[i] = std::invoke(sampler, engine);
        }
        return result;
    }

    template<typename Matrix, typename RandomEngine, typename Sampler>
        requires std::constructible_from<Matrix, std::size_t, std::size_t> &&
                 requires(Matrix &value, RandomEngine &engine, Sampler &sampler) {
                     value(0, 0) = std::invoke(sampler, engine);
                 }
    Matrix random_matrix(std::size_t rows, std::size_t columns, RandomEngine &engine, Sampler sampler) {
        Matrix result(rows, columns);
        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t column = 0; column < columns; ++column) {
                result(row, column) = std::invoke(sampler, engine);
            }
        }
        return result;
    }

    template<typename Vector, typename RandomEngine>
    Vector random_vector(std::size_t size, RandomEngine &engine) {
        using value_type = typename Vector::value_type;
        using field_type = typename value_type::field_type;
        return random_vector<Vector>(size, engine,
                                     [](auto &rng) { return value_type(algebra::random_element<field_type>(rng)); });
    }

    template<typename Matrix, typename RandomEngine>
    Matrix random_matrix(std::size_t rows, std::size_t columns, RandomEngine &engine) {
        using value_type = typename Matrix::value_type;
        using field_type = typename value_type::field_type;
        return random_matrix<Matrix>(rows, columns, engine,
                                     [](auto &rng) { return value_type(algebra::random_element<field_type>(rng)); });
    }
}    // namespace nil::crypto3::math
