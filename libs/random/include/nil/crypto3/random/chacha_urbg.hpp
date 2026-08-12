#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include <nil/crypto3/random/chacha.hpp>

namespace nil::crypto3::random {
    // Adapts the byte-oriented ChaCha generator to UniformRandomBitGenerator.
    template<typename UInt = std::uint64_t, typename ChaCha = chacha<>>
    class chacha_urbg {
        static_assert(std::is_unsigned_v<UInt>);

    public:
        using result_type = UInt;

        chacha_urbg() = default;

        template<typename SeedRange>
        explicit chacha_urbg(const SeedRange &seed) : engine_(seed) {
        }

        static constexpr result_type min() {
            return std::numeric_limits<result_type>::min();
        }

        static constexpr result_type max() {
            return std::numeric_limits<result_type>::max();
        }

        result_type operator()() {
            std::array<std::uint8_t, sizeof(result_type)> bytes {};
            engine_.generate(bytes.begin(), bytes.end());

            result_type result = 0;
            for (std::size_t i = 0; i < bytes.size(); ++i) {
                result |= static_cast<result_type>(bytes[i]) << (i * 8);
            }
            return result;
        }

    private:
        ChaCha engine_;
    };
}    // namespace nil::crypto3::random
