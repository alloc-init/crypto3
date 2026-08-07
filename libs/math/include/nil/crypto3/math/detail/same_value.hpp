#pragma once

namespace nil::crypto3::math::detail {

    template<typename T, T... v>
    struct all_same_value { };

    template<typename T, T v1, T v2, T... rest>
    struct all_same_value<T, v1, v2, rest...> : all_same_value<T, v2, rest...> {
        static_assert(v1 == v2, "All values in the template parameter list must be equal");
        static constexpr T value = v1;
    };

    template<typename T, T v>
    struct all_same_value<T, v> {
        static constexpr T value = v;
    };

    template<typename First, typename... Rest>
    constexpr bool all_same_size(const First &first, const Rest &...rest) {
        // Get the size of the first element
        constexpr const std::size_t first_size = first.size();

        // Fold expression to check all sizes are the same
        return ((rest.size() == first_size) && ...);
    }
}    // namespace nil::crypto3::math::detail
