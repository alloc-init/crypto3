//---------------------------------------------------------------------------//
// Copyright (c) 2024 Mikhail Komarov <nemo@allocin.it>
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

#ifndef CRYPTO3_RANDOM_CT_LCG_HPP
#define CRYPTO3_RANDOM_CT_LCG_HPP

#include <type_traits>

#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <nil/crypto3/algebra/type_traits.hpp>

namespace nil {
    namespace crypto3 {
        namespace random {
            template<typename Integer>
            class ct_uniform_int_distribution {
            public:
                using result_type = Integer;

                // Default constructor with range [0, 1]
                constexpr ct_uniform_int_distribution(Integer min_val = 0, Integer max_val = 1)
                        : min_(min_val), max_(max_val) {
                    if (min_ > max_) {
                        throw std::invalid_argument("Minimum value cannot be greater than maximum value.");
                    }
                }

                // Accessors for the range
                constexpr Integer min() const { return min_; }

                constexpr Integer max() const { return max_; }

                // Main distribution function to generate a number in the range [min, max]
                template<typename UniformRandomBitGenerator>
                constexpr result_type operator()(UniformRandomBitGenerator &rng) const {
                    return min_ + (rng() % (max_ - min_ + 1));
                }

            private:
                Integer min_;
                Integer max_;
            };

            template<typename Integer, Integer A, Integer C, Integer M>
            struct ct_lcg {
                using result_type = Integer;

                constexpr ct_lcg(result_type seed = 0) : state(seed) {}

                // Returns the minimum value
                static constexpr result_type min() {
                    return 0;
                }

                // Returns the maximum value
                static constexpr result_type max() {
                    return M - 1;
                }

                // Computes the nth random number given a seed
                static constexpr result_type value(result_type seed, unsigned n = 1) {
                    result_type state = seed;
                    for (unsigned i = 0; i < n; ++i) {
                        state = next(state);
                    }
                    return state;
                }

                constexpr result_type operator()() {
                    // BC: the following line causes value to run for a very very long time under certain conditions
                    state = value(state, (A * state + C) % M);

                    return state;
                }

            private:
                result_type state;

                // Computes the next state
                static constexpr result_type next(result_type current_state) {
                    return (A * current_state + C) % M;
                }
            };
        }    // namespace random
    }        // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_RANDOM_ALGEBRAIC_RANDOM_DEVICE_HPP
