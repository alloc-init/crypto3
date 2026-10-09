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

#ifndef CRYPTO3_MATH_POLYNOMIAL_BACKEND_HPP
#define CRYPTO3_MATH_POLYNOMIAL_BACKEND_HPP

#include <concepts>
#include <cstddef>
#include <limits>
#include <utility>

#include <nil/crypto3/math/polynomial/concepts.hpp>
#include <nil/crypto3/math/polynomial/polynomial_observation.hpp>

namespace nil::crypto3::math::polynomial_arithmetic {

    /** Algorithm-selection parameters shared by higher-level polynomial operations. */
    struct polynomial_context_options {
        // Long division is selected when either enabled coefficient-count cutoff matches.
        std::size_t basecase_divisor_coefficient_cutoff = 10;
        std::size_t basecase_quotient_coefficient_cutoff = 2;

        // Recursive half-GCD switches to iterative matrix construction at or below this coefficient count. Zero
        // disables the iterative base case.
        std::size_t half_gcd_basecase_cutoff = 30;

        // GCD uses half-GCD when its smaller operand has at least this many coefficients. Zero disables half-GCD;
        // callers should select this backend-dependent crossover from benchmarks.
        std::size_t gcd_half_gcd_cutoff = 0;

        // Brent-Kung modular composition caches at most this many reduced powers of the inner polynomial. Each power
        // has at most degree(divisor) coefficients. The limit must be positive; the default permits the algorithm's
        // square-root block size.
        std::size_t modular_composition_cached_power_limit = std::numeric_limits<std::size_t>::max();
    };

    /**
     * Interface for interchangeable polynomial multiplication implementations.
     * Higher-level polynomial algorithms use these operations without depending
     * on how products are computed.
     *
     * The associated polynomial type must use the coefficient representation.
     * Inputs use canonical form and every operation produces canonical
     * output. The output may alias either input. multiply_low computes the product
     * modulo X^coefficient_count; when coefficient_count is zero, it stores [0].
     *
     * Operations are invoked on a mutable backend so implementations may update
     * reusable caches or scratch storage; stateless implementations may still
     * declare their operations const.
     */
    template<typename Backend>
    concept PolynomialBackend =
        requires {
            typename Backend::polynomial_type;
            requires CoefficientPolynomial<typename Backend::polynomial_type>;
        } && requires(Backend &backend, typename Backend::polynomial_type &output,
                      const typename Backend::polynomial_type &left, const typename Backend::polynomial_type &right,
                      std::size_t coefficient_count) {
            { backend.multiply(output, left, right) } -> std::same_as<void>;
            { backend.square(output, left) } -> std::same_as<void>;
            { backend.multiply_low(output, left, right, coefficient_count) } -> std::same_as<void>;
        };

    /**
     * Owns the multiplication implementation and algorithm-selection parameters
     * used by a sequence of higher-level polynomial operations. Keeping one backend
     * alive lets those operations reuse implementation-specific configuration,
     * precomputed state, and scratch storage instead of rebuilding them for every
     * product.
     *
     * Higher-level algorithms use the context without managing backend-specific
     * plans, configuration, or scratch storage. A context can be reused sequentially,
     * but callers must use separate contexts for concurrent operations.
     *
     * An enabled observer is borrowed and must outlive the context and its scopes. Copies of a context borrow
     * the same observer; use separate observers for concurrent workers. Default observation has empty storage
     * and requires no observer methods on the backend. It does not change optional backend capabilities.
     * The context emits operation scopes. A backend may additionally accept the borrowed observer as its last
     * argument to emit kernel phases; otherwise its original operation runs inside the context's scope. These
     * extra overloads are selected only when kernel observation is enabled and must preserve the original result,
     * acceptance, and exception semantics. Kernel-only observers receive no algorithm scopes.
     */
    template<PolynomialBackend Backend, PolynomialObserver Observer = no_polynomial_observer>
    class polynomial_context {
    public:
        using backend_type = Backend;
        using polynomial_type = typename backend_type::polynomial_type;
        using value_type = typename polynomial_type::value_type;
        using options_type = polynomial_context_options;
        using observer_type = Observer;

        polynomial_context()
            requires(!(Observer::observe_algorithms || Observer::observe_kernels))
        = default;

        explicit polynomial_context(backend_type backend, options_type options = {})
            requires(!(Observer::observe_algorithms || Observer::observe_kernels))
            : backend_(std::move(backend)), options_(options) {
        }

        polynomial_context(backend_type backend, observer_type &observer, options_type options = {}) :
            backend_(std::move(backend)), options_(options), observer_(observer) {
        }

        // Metadata is built only for enabled stages. Put any observation-only computation in this factory.
        template<polynomial_stage Stage, detail::PolynomialMetadataFactory Factory>
        [[nodiscard]] auto observe(Factory &&metadata) noexcept {
            using scope_type = polynomial_observation_scope<observer_type, Stage>;
            if constexpr (scope_type::enabled) {
                return observe_polynomial<Stage>(*observer_.observer, std::forward<Factory>(metadata));
            } else {
                return scope_type();
            }
        }

        template<polynomial_stage Stage>
        [[nodiscard]] auto observe() noexcept {
            return observe<Stage>([]() noexcept { return polynomial_metadata {}; });
        }

        void multiply(polynomial_type &output, const polynomial_type &left, const polynomial_type &right) {
            auto scope = observe<polynomial_stage::multiply>([&]() noexcept {
                return polynomial_metadata {{{polynomial_metric::input_coefficients, left.size()},
                                             {polynomial_metric::second_input_coefficients, right.size()}}};
            });
            if constexpr (Observer::observe_kernels && requires(observer_type &observer) {
                              { backend_.multiply(output, left, right, observer) } -> std::same_as<void>;
                          }) {
                backend_.multiply(output, left, right, *observer_.observer);
                scope.set_result_metadata([&]() noexcept {
                    return polynomial_metadata {{{polynomial_metric::output_coefficients, output.size()},
                                                 {polynomial_metric::detailed_backend, 1}}};
                });
            } else {
                backend_.multiply(output, left, right);
                scope.set_result_metadata([&]() noexcept {
                    return polynomial_metadata {{{polynomial_metric::output_coefficients, output.size()},
                                                 {polynomial_metric::detailed_backend, 0}}};
                });
            }
        }

        void square(polynomial_type &output, const polynomial_type &input) {
            auto scope = observe<polynomial_stage::square>([&]() noexcept {
                return polynomial_metadata {{{polynomial_metric::input_coefficients, input.size()}}};
            });
            if constexpr (Observer::observe_kernels && requires(observer_type &observer) {
                              { backend_.square(output, input, observer) } -> std::same_as<void>;
                          }) {
                backend_.square(output, input, *observer_.observer);
                scope.set_result_metadata([&]() noexcept {
                    return polynomial_metadata {{{polynomial_metric::output_coefficients, output.size()},
                                                 {polynomial_metric::detailed_backend, 1}}};
                });
            } else {
                backend_.square(output, input);
                scope.set_result_metadata([&]() noexcept {
                    return polynomial_metadata {{{polynomial_metric::output_coefficients, output.size()},
                                                 {polynomial_metric::detailed_backend, 0}}};
                });
            }
        }

        void multiply_low(polynomial_type &output, const polynomial_type &left, const polynomial_type &right,
                          std::size_t coefficient_count) {
            auto scope = observe<polynomial_stage::multiply_low>([&]() noexcept {
                return polynomial_metadata {{{polynomial_metric::input_coefficients, left.size()},
                                             {polynomial_metric::second_input_coefficients, right.size()},
                                             {polynomial_metric::low_coefficient_count, coefficient_count}}};
            });
            if constexpr (Observer::observe_kernels && requires(observer_type &observer) {
                              {
                                  backend_.multiply_low(output, left, right, coefficient_count, observer)
                              } -> std::same_as<void>;
                          }) {
                backend_.multiply_low(output, left, right, coefficient_count, *observer_.observer);
                scope.set_result_metadata([&]() noexcept {
                    return polynomial_metadata {{{polynomial_metric::output_coefficients, output.size()},
                                                 {polynomial_metric::detailed_backend, 1}}};
                });
            } else {
                backend_.multiply_low(output, left, right, coefficient_count);
                scope.set_result_metadata([&]() noexcept {
                    return polynomial_metadata {{{polynomial_metric::output_coefficients, output.size()},
                                                 {polynomial_metric::detailed_backend, 0}}};
                });
            }
        }

        // Optional capabilities: backends without prepared low products still satisfy PolynomialBackend.
        // Forward the backend's owned preparation without adding another context or exposing its scratch storage.
        auto prepare_low_product(const polynomial_type &fixed, std::size_t variable_coefficient_count,
                                 std::size_t coefficient_count)
            requires requires(backend_type &backend, const polynomial_type &operand, std::size_t count) {
                backend.prepare_low_product(operand, count, count);
            }
        {
            if constexpr (!Observer::observe_kernels) {
                return backend_.prepare_low_product(fixed, variable_coefficient_count, coefficient_count);
            } else {
                auto scope = observe<polynomial_stage::prepare_low_product>([&]() noexcept {
                    return polynomial_metadata {
                        {{polynomial_metric::input_coefficients, fixed.size()},
                         {polynomial_metric::second_input_coefficients, variable_coefficient_count},
                         {polynomial_metric::low_coefficient_count, coefficient_count}}};
                });
                auto prepared = [&]() {
                    if constexpr (requires(observer_type &observer) {
                                      backend_.prepare_low_product(fixed, variable_coefficient_count, coefficient_count,
                                                                   observer);
                                  }) {
                        return backend_.prepare_low_product(fixed, variable_coefficient_count, coefficient_count,
                                                            *observer_.observer);
                    } else {
                        return backend_.prepare_low_product(fixed, variable_coefficient_count, coefficient_count);
                    }
                }();
                if (!prepared) {
                    scope.decline();
                }
                scope.set_result_metadata([&]() noexcept {
                    polynomial_metadata result {{{polynomial_metric::accepted, bool(prepared)}}};
                    if (prepared) {
                        if constexpr (requires {
                                          { prepared->transform_size() } noexcept -> std::convertible_to<std::size_t>;
                                      }) {
                            result[1] = {polynomial_metric::transform_length, prepared->transform_size()};
                        }
                        if constexpr (requires {
                                          { prepared->storage_bytes() } noexcept -> std::convertible_to<std::size_t>;
                                      }) {
                            result[2] = {polynomial_metric::prepared_storage_bytes, prepared->storage_bytes()};
                        }
                    }
                    return result;
                });
                return prepared;
            }
        }

        template<typename PreparedOperand>
        bool try_multiply_low_prepared(polynomial_type &output, const polynomial_type &left,
                                       const PreparedOperand &prepared, std::size_t coefficient_count)
            requires requires(backend_type &backend, polynomial_type &result, const polynomial_type &operand,
                              const PreparedOperand &fixed, std::size_t count) {
                { backend.try_multiply_low_prepared(result, operand, fixed, count) } -> std::same_as<bool>;
            }
        {
            if constexpr (!Observer::observe_kernels) {
                return backend_.try_multiply_low_prepared(output, left, prepared, coefficient_count);
            } else {
                auto scope = observe<polynomial_stage::try_multiply_low_prepared>([&]() noexcept {
                    return polynomial_metadata {{{polynomial_metric::input_coefficients, left.size()},
                                                 {polynomial_metric::low_coefficient_count, coefficient_count}}};
                });
                const bool accepted = [&]() {
                    if constexpr (requires(observer_type &observer) {
                                      {
                                          backend_.try_multiply_low_prepared(output, left, prepared, coefficient_count,
                                                                             observer)
                                      } -> std::same_as<bool>;
                                  }) {
                        return backend_.try_multiply_low_prepared(output, left, prepared, coefficient_count,
                                                                  *observer_.observer);
                    } else {
                        return backend_.try_multiply_low_prepared(output, left, prepared, coefficient_count);
                    }
                }();
                if (!accepted) {
                    scope.decline();
                }
                scope.set_result_metadata([&]() noexcept {
                    polynomial_metadata result {{{polynomial_metric::accepted, accepted}}};
                    if (accepted) {
                        result[1] = {polynomial_metric::output_coefficients, output.size()};
                    }
                    return result;
                });
                return accepted;
            }
        }

        auto prepare_cyclic_remainder(const polynomial_type &divisor, std::size_t quotient_coefficient_count)
            requires requires(backend_type &backend, const polynomial_type &operand, std::size_t count) {
                backend.prepare_cyclic_remainder(operand, count);
            }
        {
            if constexpr (!Observer::observe_kernels) {
                return backend_.prepare_cyclic_remainder(divisor, quotient_coefficient_count);
            } else {
                auto scope = observe<polynomial_stage::prepare_cyclic_remainder>([&]() noexcept {
                    return polynomial_metadata {
                        {{polynomial_metric::input_coefficients, divisor.size()},
                         {polynomial_metric::quotient_coefficients, quotient_coefficient_count}}};
                });
                auto prepared = [&]() {
                    if constexpr (requires(observer_type &observer) {
                                      backend_.prepare_cyclic_remainder(divisor, quotient_coefficient_count, observer);
                                  }) {
                        return backend_.prepare_cyclic_remainder(divisor, quotient_coefficient_count,
                                                                 *observer_.observer);
                    } else {
                        return backend_.prepare_cyclic_remainder(divisor, quotient_coefficient_count);
                    }
                }();
                if (!prepared) {
                    scope.decline();
                }
                scope.set_result_metadata([&]() noexcept {
                    polynomial_metadata result {{{polynomial_metric::accepted, bool(prepared)}}};
                    if (prepared) {
                        if constexpr (requires {
                                          { prepared->transform_size() } noexcept -> std::convertible_to<std::size_t>;
                                      }) {
                            result[1] = {polynomial_metric::transform_length, prepared->transform_size()};
                        }
                        if constexpr (requires {
                                          { prepared->storage_bytes() } noexcept -> std::convertible_to<std::size_t>;
                                      }) {
                            result[2] = {polynomial_metric::prepared_storage_bytes, prepared->storage_bytes()};
                        }
                    }
                    return result;
                });
                return prepared;
            }
        }

        template<typename PreparedRemainder>
        bool try_cyclic_remainder(polynomial_type &output, const polynomial_type &dividend,
                                  const polynomial_type &quotient, const PreparedRemainder &prepared)
            requires requires(backend_type &backend, polynomial_type &result, const polynomial_type &operand,
                              const PreparedRemainder &fixed) {
                { backend.try_cyclic_remainder(result, operand, operand, fixed) } -> std::same_as<bool>;
            }
        {
            if constexpr (!Observer::observe_kernels) {
                return backend_.try_cyclic_remainder(output, dividend, quotient, prepared);
            } else {
                auto scope = observe<polynomial_stage::try_cyclic_remainder>([&]() noexcept {
                    return polynomial_metadata {{{polynomial_metric::input_coefficients, dividend.size()},
                                                 {polynomial_metric::quotient_coefficients, quotient.size()}}};
                });
                const bool accepted = [&]() {
                    if constexpr (requires(observer_type &observer) {
                                      {
                                          backend_.try_cyclic_remainder(output, dividend, quotient, prepared, observer)
                                      } -> std::same_as<bool>;
                                  }) {
                        return backend_.try_cyclic_remainder(output, dividend, quotient, prepared, *observer_.observer);
                    } else {
                        return backend_.try_cyclic_remainder(output, dividend, quotient, prepared);
                    }
                }();
                if (!accepted) {
                    scope.decline();
                }
                scope.set_result_metadata([&]() noexcept {
                    polynomial_metadata result {{{polynomial_metric::accepted, accepted}}};
                    if (accepted) {
                        result[1] = {polynomial_metric::output_coefficients, output.size()};
                    }
                    return result;
                });
                return accepted;
            }
        }

        const options_type &options() const {
            return options_;
        }

    private:
        backend_type backend_;
        options_type options_;
        [[no_unique_address]] detail::polynomial_observer_reference<observer_type> observer_;
    };

}    // namespace nil::crypto3::math::polynomial_arithmetic

#endif    // CRYPTO3_MATH_POLYNOMIAL_BACKEND_HPP
