//---------------------------------------------------------------------------//
// Copyright (c) 2026 Alloc Init
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

#define BOOST_TEST_MODULE modified_sqap_prover_test

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/prover.hpp>

#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/operations/lagrange_interpolation.hpp>
#include <nil/crypto3/random/chacha_urbg.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/generator.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/policy.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sqap/transcript.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_value_type = curve_type::scalar_field_type::value_type;
    using base_value_type = curve_type::base_field_type::value_type;
    using g1_value_type = curve_type::g1_type<>::value_type;
    using pairing_policy_type = nil::crypto3::zk::snark::modified_sqap_bn254_exact_pairing_policy;
    using transcript_policy_type = nil::crypto3::zk::snark::modified_sqap_bn254_poseidon_transcript_policy;
    using policy_type =
        nil::crypto3::zk::snark::modified_sqap_policy<curve_type, pairing_policy_type, transcript_policy_type>;
    using prover_type = nil::crypto3::zk::snark::modified_sqap_prover<policy_type>;
    using generator_type = nil::crypto3::zk::snark::modified_sqap_generator<policy_type>;
    using deterministic_generator_type =
        nil::crypto3::zk::snark::detail::modified_sqap_deterministic_generator<policy_type>;
    using reduction_type =
        nil::crypto3::zk::snark::reductions::modified_sqap_to_polynomials<curve_type::scalar_field_type>;
    using polynomial_type = reduction_type::polynomial_type;
    using auxiliary_input_type = policy_type::auxiliary_input_type;
    using proving_key_type = policy_type::proving_key_type;
    using system_type = policy_type::constraint_system_type;
    using constraint_type = system_type::constraint_type;
    using variable_type = constraint_type::variable_type;
    using combination_type = constraint_type::linear_combination_type;

    combination_type raw_combination(std::initializer_list<std::pair<std::size_t, int>> terms) {
        combination_type result;
        for (const auto &[index, coefficient] : terms) {
            result.add_term(variable_type(index), scalar_value_type(coefficient));
        }
        return result;
    }

    constraint_type binding_row() {
        return {{}, combination_type(-variable_type(0))};
    }

    system_type test_system() {
        return {3,
                {binding_row(),
                 {raw_combination({{1, 2}, {1, -1}}), raw_combination({{2, 1}})},
                 {raw_combination({{2, 1}}), raw_combination({{2, -2}})}}};
    }

    proving_key_type generate_proving_key(const system_type &source) {
        // Fixed seed for reproducible tests with the production setup and transcript.
        const std::array<std::uint8_t, 32> seed = {};
        nil::crypto3::random::chacha_urbg<> random_source(seed);
        return generator_type::process(source, random_source).first;
    }

    deterministic_generator_type::trapdoor_type test_trapdoor() {
        return {scalar_value_type(2),
                scalar_value_type(3),
                {scalar_value_type(5), scalar_value_type(7), scalar_value_type(11), scalar_value_type(13)}};
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sqap_prover_test_suite)

BOOST_AUTO_TEST_CASE(generated_keys_validate_without_modification) {
    for (const std::size_t logical_rows : {1, 2, 3, 4, 5}) {
        BOOST_TEST_CONTEXT("logical rows: " << logical_rows) {
            auto source = test_system();
            source.constraints.resize(logical_rows, binding_row());
            const auto key = generate_proving_key(source);
            const auto original = key;

            BOOST_CHECK_NO_THROW(prover_type::validate(key));
            BOOST_CHECK(key == original);
        }
    }
}

BOOST_AUTO_TEST_CASE(minimum_domain_accepts_only_an_empty_h_opening_query) {
    auto key = generate_proving_key({1, {binding_row()}});
    BOOST_REQUIRE_EQUAL(key.verification_key.domain_size, 2);
    BOOST_REQUIRE(key.T_H.empty());
    BOOST_CHECK_NO_THROW(prover_type::validate(key));

    key.T_H.push_back(g1_value_type::zero());
    BOOST_CHECK_THROW(prover_type::validate(key), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(rejects_short_and_long_queries) {
    const auto original = generate_proving_key(test_system());
    const std::array queries = {
        std::pair {"W", &proving_key_type::W},     std::pair {"H_query", &proving_key_type::H_query},
        std::pair {"T_A", &proving_key_type::T_A}, std::pair {"T_C", &proving_key_type::T_C},
        std::pair {"T_H", &proving_key_type::T_H}, std::pair {"T_Z", &proving_key_type::T_Z}};

    for (const auto &[name, member] : queries) {
        BOOST_TEST_CONTEXT("query: " << name) {
            BOOST_REQUIRE(!(original.*member).empty());
            auto shorter = original;
            (shorter.*member).pop_back();
            BOOST_CHECK_THROW(prover_type::validate(shorter), std::invalid_argument);

            auto longer = original;
            (longer.*member).push_back(g1_value_type::zero());
            BOOST_CHECK_THROW(prover_type::validate(longer), std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_inconsistent_verification_key_dimensions) {
    const auto original = generate_proving_key(test_system());
    for (const auto n : {std::size_t(0), std::size_t(2), std::size_t(4), std::numeric_limits<std::size_t>::max()}) {
        BOOST_TEST_CONTEXT("num_variables: " << n) {
            auto key = original;
            key.verification_key.num_variables = n;
            BOOST_CHECK_THROW(prover_type::validate(key), std::invalid_argument);
        }
    }
    for (const auto m : {std::size_t(0), std::size_t(1), std::size_t(2), std::size_t(3), std::size_t(5), std::size_t(8),
                         std::numeric_limits<std::size_t>::max()}) {
        BOOST_TEST_CONTEXT("domain_size: " << m) {
            auto key = original;
            key.verification_key.domain_size = m;
            BOOST_CHECK_THROW(prover_type::validate(key), std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_invalid_constraint_systems) {
    BOOST_CHECK_THROW(prover_type::validate(proving_key_type {}), std::invalid_argument);
    const auto original = generate_proving_key(test_system());

    auto zero_variables = original;
    zero_variables.constraint_system.witness_size = 0;
    BOOST_CHECK_THROW(prover_type::validate(zero_variables), std::invalid_argument);

    auto empty_rows = original;
    empty_rows.constraint_system.constraints.clear();
    BOOST_CHECK_THROW(prover_type::validate(empty_rows), std::invalid_argument);

    auto missing_binding = original;
    missing_binding.constraint_system.constraints.front() = original.constraint_system.constraints[1];
    BOOST_CHECK_THROW(prover_type::validate(missing_binding), std::invalid_argument);

    auto out_of_range = original;
    out_of_range.constraint_system.constraints[1].a.add_term(
        variable_type(out_of_range.constraint_system.num_variables()), scalar_value_type::zero());
    BOOST_CHECK_THROW(prover_type::validate(out_of_range), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(rejects_changed_circuit_dimensions) {
    const auto original = generate_proving_key(test_system());

    auto more_variables = original;
    ++more_variables.constraint_system.witness_size;
    BOOST_CHECK_THROW(prover_type::validate(more_variables), std::invalid_argument);

    auto more_rows = original;
    more_rows.constraint_system.constraints.resize(5, binding_row());
    BOOST_CHECK_THROW(prover_type::validate(more_rows), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(rejects_changed_digest_or_same_shape_circuit) {
    const auto original = generate_proving_key(test_system());

    auto changed_digest = original;
    changed_digest.verification_key.circuit_digest += base_value_type::one();
    BOOST_CHECK_THROW(prover_type::validate(changed_digest), std::invalid_argument);

    auto changed_coefficient = original;
    changed_coefficient.constraint_system.constraints[1].a.terms.front().coeff += scalar_value_type::one();
    BOOST_REQUIRE(changed_coefficient.constraint_system.is_valid());
    BOOST_CHECK_THROW(prover_type::validate(changed_coefficient), std::invalid_argument);

    auto changed_index = original;
    changed_index.constraint_system.constraints[1].a.terms.front().index = 2;
    BOOST_REQUIRE(changed_index.constraint_system.is_valid());
    BOOST_CHECK_THROW(prover_type::validate(changed_index), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(digest_binds_logical_rows_and_witness_dimension) {
    const auto original = generate_proving_key(test_system());

    // A fourth logical row preserves the padded domain size and every query length.
    auto more_rows = original;
    more_rows.constraint_system.constraints.push_back(binding_row());
    BOOST_REQUIRE(more_rows.constraint_system.is_valid());
    BOOST_CHECK_THROW(prover_type::validate(more_rows), std::invalid_argument);

    // Matching the new witness dimension in the metadata and W leaves the digest stale.
    auto more_variables = original;
    ++more_variables.constraint_system.witness_size;
    ++more_variables.verification_key.num_variables;
    more_variables.W.push_back(g1_value_type::zero());
    BOOST_REQUIRE(more_variables.constraint_system.is_valid());
    BOOST_CHECK_THROW(prover_type::validate(more_variables), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(commitment_matches_independent_scalar_reference) {
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const auto original_witness = witness;
    const auto trapdoor = test_trapdoor();

    for (const std::size_t logical_rows : {1, 2, 3, 4, 5}) {
        BOOST_TEST_CONTEXT("logical rows: " << logical_rows) {
            auto source = test_system();
            source.constraints.resize(logical_rows, binding_row());
            const auto key =
                deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor)
                    .first;
            const auto original_key = key;
            const auto result = prover_type::commit(key, u, witness);

            const auto domain = reduction_type::get_domain(source);
            const auto m = domain->size();
            BOOST_REQUIRE_EQUAL(result.polynomials.A.size(), m);
            BOOST_REQUIRE_EQUAL(result.polynomials.C.size(), m);
            BOOST_REQUIRE_EQUAL(result.polynomials.H.size(), m - 1);

            // Reference interpolation uses products of linear factors instead of inverse FFTs.
            std::vector<std::pair<scalar_value_type, scalar_value_type>> a_points, c_points;
            for (std::size_t j = 0; j < m; ++j) {
                const auto x = domain->get_domain_element(j);
                const auto a = j < logical_rows ? source.constraints[j].a.evaluate(witness) : scalar_value_type::zero();
                const auto c = j < logical_rows ? source.constraints[j].c.evaluate(witness) : -u;
                a_points.emplace_back(x, a);
                c_points.emplace_back(x, c);
            }
            auto expected_a = nil::crypto3::math::lagrange_interpolation(a_points);
            auto expected_c = nil::crypto3::math::lagrange_interpolation(c_points);
            expected_a.resize(m, scalar_value_type::zero());
            expected_c.resize(m, scalar_value_type::zero());
            BOOST_CHECK(result.polynomials.A == expected_a);
            BOOST_CHECK(result.polynomials.C == expected_c);

            polynomial_type expected_z(m + 1, scalar_value_type::zero());
            expected_z[0] = -scalar_value_type::one();
            expected_z[m] = scalar_value_type::one();
            BOOST_CHECK(result.Z == expected_z);

            // Recover H(tau) directly from the relation, without using the computed H coefficients or queries.
            const auto a_at_tau = expected_a.evaluate(trapdoor.tau);
            const auto c_at_tau = expected_c.evaluate(trapdoor.tau);
            const auto h_at_tau =
                (a_at_tau.squared() - c_at_tau - u) * (trapdoor.tau.pow(m) - scalar_value_type::one()).inversed();
            BOOST_CHECK_EQUAL(result.polynomials.H.evaluate(trapdoor.tau), h_at_tau);
            const auto expected_scalar = trapdoor.gamma * (trapdoor.alpha[0] * a_at_tau + trapdoor.alpha[1] * c_at_tau +
                                                           trapdoor.alpha[2] * h_at_tau);
            BOOST_CHECK_EQUAL(result.P, expected_scalar * g1_value_type::one());
            BOOST_CHECK(key == original_key);
            BOOST_CHECK(witness == original_witness);
        }
    }
}

BOOST_AUTO_TEST_CASE(commitment_handles_zero_witness_entries_and_zero_h) {
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type::zero()};
    const auto trapdoor = test_trapdoor();
    for (const std::size_t logical_rows : {2, 3, 5}) {
        BOOST_TEST_CONTEXT("logical rows: " << logical_rows) {
            system_type source {2, {binding_row(), {raw_combination({{1, 1}}), raw_combination({{0, -1}})}}};
            source.constraints.resize(logical_rows, binding_row());
            const auto key =
                deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor)
                    .first;
            BOOST_REQUIRE(!key.W[1].is_zero());
            BOOST_REQUIRE(!key.H_query.front().is_zero());

            const auto result = prover_type::commit(key, u, witness);
            BOOST_CHECK(result.polynomials.A.is_zero());
            BOOST_CHECK(result.polynomials.H.is_zero());
            // Every row evaluates to A = 0 and C = -u, so H = 0 and P = [-gamma * alpha_C * u]_1.
            const auto expected_scalar = -trapdoor.gamma * trapdoor.alpha[1] * u;
            BOOST_CHECK_EQUAL(result.P, expected_scalar * g1_value_type::one());
        }
    }
}

BOOST_AUTO_TEST_CASE(commitment_accepts_identity_query_points) {
    const auto source = test_system();
    auto trapdoor = test_trapdoor();
    trapdoor.alpha[0] = scalar_value_type::zero();
    trapdoor.alpha[1] = scalar_value_type::zero();
    trapdoor.alpha[2] = scalar_value_type::zero();
    const auto key =
        deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor).first;
    for (const auto &point : key.W) {
        BOOST_REQUIRE(point.is_zero());
    }
    for (const auto &point : key.H_query) {
        BOOST_REQUIRE(point.is_zero());
    }

    const scalar_value_type u(3);
    const auto result = prover_type::commit(key, u, {u, scalar_value_type(2), scalar_value_type(1)});
    BOOST_REQUIRE(!result.polynomials.H.is_zero());
    BOOST_CHECK(result.P.is_zero());
}

BOOST_AUTO_TEST_CASE(commitment_allows_cancellation_to_identity) {
    auto source = test_system();
    source.constraints.resize(2);
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};

    // A = 1 - X, C = -1 - 2X and H = 1. At tau = 2 the commitment scalar is
    // gamma * (-alpha_A - 5 * alpha_C + alpha_H).
    const std::array weights = {
        std::pair {scalar_value_type(7), scalar_value_type(40)},     // The W and H sums cancel.
        std::pair {-scalar_value_type(1), scalar_value_type(0)}};    // The W sum cancels internally.
    for (const auto &[alpha_c, alpha_h] : weights) {
        BOOST_TEST_CONTEXT("alpha_C: " << alpha_c << ", alpha_H: " << alpha_h) {
            auto trapdoor = test_trapdoor();
            trapdoor.alpha[1] = alpha_c;
            trapdoor.alpha[2] = alpha_h;
            const auto key =
                deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor)
                    .first;
            for (const auto &point : key.W) {
                BOOST_REQUIRE(!point.is_zero());
            }

            const auto result = prover_type::commit(key, u, witness);
            BOOST_REQUIRE_EQUAL(result.polynomials.H.size(), 1);
            BOOST_CHECK_EQUAL(result.polynomials.H[0], scalar_value_type::one());
            BOOST_CHECK(result.P.is_zero());
        }
    }
}

BOOST_AUTO_TEST_CASE(commitment_rejects_inconsistent_keys) {
    const auto original = generate_proving_key(test_system());
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    for (const auto member : {&proving_key_type::W, &proving_key_type::H_query}) {
        auto key = original;
        (key.*member).pop_back();
        BOOST_CHECK_THROW(prover_type::commit(key, u, witness), std::invalid_argument);
    }

    auto changed_digest = original;
    changed_digest.verification_key.circuit_digest += base_value_type::one();
    BOOST_CHECK_THROW(prover_type::commit(changed_digest, u, witness), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(commitment_rejects_invalid_public_inputs_and_witnesses) {
    const auto key = generate_proving_key(test_system());
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};

    for (const std::size_t size : {0, 2, 4}) {
        auto wrong_length = witness;
        wrong_length.resize(size, scalar_value_type::zero());
        BOOST_CHECK_THROW(prover_type::commit(key, u, wrong_length), std::invalid_argument);
    }
    auto wrong_binding = witness;
    wrong_binding[0] = scalar_value_type(5);
    BOOST_CHECK_THROW(prover_type::commit(key, u, wrong_binding), std::invalid_argument);
    BOOST_CHECK_THROW(prover_type::commit(key, scalar_value_type(5), witness), std::invalid_argument);

    auto unsatisfied = witness;
    unsatisfied[1] += scalar_value_type::one();
    BOOST_CHECK_THROW(prover_type::commit(key, u, unsatisfied), std::invalid_argument);

    // Binding-only witnesses satisfy the rows for any u, isolating the canonical parity check.
    const auto binding_key = generate_proving_key({1, {binding_row()}});
    for (const auto &even_u : {scalar_value_type::zero(), scalar_value_type(2), scalar_value_type(4)}) {
        BOOST_CHECK_THROW(prover_type::commit(binding_key, even_u, {even_u}), std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(opening_matches_explicit_cubic_example) {
    // U(X) = 5 + 2X + 3X^2 + 4X^3.
    const polynomial_type polynomial {scalar_value_type(5), scalar_value_type(2), scalar_value_type(3),
                                      scalar_value_type(4)};
    const auto original = polynomial;
    for (const auto &z :
         {scalar_value_type::zero(), scalar_value_type::one(), scalar_value_type(2), -scalar_value_type(2)}) {
        BOOST_TEST_CONTEXT("challenge: " << z) {
            const auto result = prover_type::open_polynomial(polynomial, z);
            const auto expected_evaluation = scalar_value_type(5) + scalar_value_type(2) * z +
                                             scalar_value_type(3) * z.squared() +
                                             scalar_value_type(4) * z.squared() * z;
            const polynomial_type expected_quotient {
                scalar_value_type(2) + scalar_value_type(3) * z + scalar_value_type(4) * z.squared(),
                scalar_value_type(3) + scalar_value_type(4) * z, scalar_value_type(4)};
            BOOST_CHECK_EQUAL(result.evaluation, expected_evaluation);
            BOOST_CHECK(result.quotient == expected_quotient);
            BOOST_CHECK(polynomial == original);
        }
    }
}

BOOST_AUTO_TEST_CASE(opening_of_constant_and_zero_polynomials_has_zero_quotient) {
    for (const auto &constant : {scalar_value_type::zero(), scalar_value_type(7)}) {
        for (const std::size_t coefficient_count : {1, 5}) {
            polynomial_type polynomial(coefficient_count, scalar_value_type::zero());
            polynomial[0] = constant;
            const auto original = polynomial;
            for (const auto &z : {scalar_value_type::zero(), scalar_value_type::one(), scalar_value_type(2)}) {
                BOOST_TEST_CONTEXT("constant: " << constant << ", coefficient count: " << coefficient_count
                                                << ", challenge: " << z) {
                    const auto result = prover_type::open_polynomial(polynomial, z);
                    BOOST_CHECK_EQUAL(result.evaluation, constant);
                    BOOST_REQUIRE_EQUAL(result.quotient.size(), 1);
                    BOOST_CHECK(result.quotient[0].is_zero());
                    BOOST_CHECK(polynomial == original);
                }
            }
        }
    }

    polynomial_type empty;
    empty.clear();
    const auto result = prover_type::open_polynomial(empty, scalar_value_type::zero());
    BOOST_CHECK(result.evaluation.is_zero());
    BOOST_REQUIRE_EQUAL(result.quotient.size(), 1);
    BOOST_CHECK(result.quotient[0].is_zero());
    BOOST_CHECK(empty.empty());
}

BOOST_AUTO_TEST_CASE(opening_preserves_trailing_zero_storage) {
    polynomial_type polynomial {scalar_value_type(5), scalar_value_type(2), scalar_value_type(3), scalar_value_type(4)};
    polynomial.resize(8, scalar_value_type::zero());
    const auto original = polynomial;

    const auto result = prover_type::open_polynomial(polynomial, scalar_value_type(2));
    const polynomial_type expected_quotient {scalar_value_type(24), scalar_value_type(11), scalar_value_type(4)};
    BOOST_CHECK_EQUAL(result.evaluation, scalar_value_type(53));
    BOOST_CHECK(result.quotient == expected_quotient);
    BOOST_CHECK_EQUAL(polynomial.size(), 8);
    BOOST_CHECK(polynomial == original);
}

BOOST_AUTO_TEST_CASE(openings_satisfy_coefficient_identities_at_zero_and_domain_points) {
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const nil::crypto3::math::polynomial_arithmetic::schoolbook_backend<scalar_value_type> reference;

    for (const std::size_t logical_rows : {1, 2, 3, 5}) {
        auto source = test_system();
        source.constraints.resize(logical_rows, binding_row());
        const auto key = generate_proving_key(source);
        const auto committed = prover_type::commit(key, u, witness);
        const auto domain = reduction_type::get_domain(source);
        const auto m = domain->size();
        std::vector<scalar_value_type> challenges {scalar_value_type::zero(), scalar_value_type(2)};
        for (std::size_t j = 0; j < m; ++j) {
            challenges.push_back(domain->get_domain_element(j));
        }
        const std::array polynomials = {&committed.polynomials.A, &committed.polynomials.C, &committed.polynomials.H,
                                        &committed.Z};
        const std::array<const char *, 4> names = {"A", "C", "H", "Z"};
        const std::array<std::size_t, 4> quotient_lengths = {m - 1, m - 1, m - 2, m};

        for (const auto &z : challenges) {
            for (std::size_t i = 0; i < polynomials.size(); ++i) {
                BOOST_TEST_CONTEXT("logical rows: " << logical_rows << ", polynomial: " << names[i]
                                                    << ", challenge: " << z) {
                    const auto original = *polynomials[i];
                    const auto result = prover_type::open_polynomial(*polynomials[i], z);
                    if (result.quotient.is_zero()) {
                        BOOST_CHECK_EQUAL(result.quotient.size(), 1);
                    } else {
                        BOOST_CHECK_LE(result.quotient.size(), quotient_lengths[i]);
                        BOOST_CHECK(!result.quotient.back().is_zero());
                    }

                    // Independently reconstruct U = (X - z) * quotient + evaluation by coefficient multiplication.
                    const polynomial_type divisor {-z, scalar_value_type::one()};
                    polynomial_type reconstructed;
                    reference.multiply(reconstructed, divisor, result.quotient);
                    reconstructed[0] += result.evaluation;
                    reconstructed.condense();
                    auto expected = original;
                    expected.condense();
                    BOOST_CHECK(reconstructed == expected);
                    BOOST_CHECK(*polynomials[i] == original);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(minimum_domain_h_opening_has_zero_quotient) {
    auto source = test_system();
    source.constraints.resize(2);
    const auto key = generate_proving_key(source);
    const scalar_value_type u(3);
    const auto committed = prover_type::commit(key, u, {u, scalar_value_type(2), scalar_value_type(1)});
    BOOST_REQUIRE_EQUAL(key.verification_key.domain_size, 2);
    BOOST_REQUIRE(key.T_H.empty());
    BOOST_REQUIRE_EQUAL(committed.polynomials.H.size(), 1);
    BOOST_REQUIRE_EQUAL(committed.polynomials.H[0], scalar_value_type::one());

    for (const auto &z :
         {scalar_value_type::zero(), scalar_value_type::one(), -scalar_value_type::one(), scalar_value_type(2)}) {
        const auto result = prover_type::open_polynomial(committed.polynomials.H, z);
        BOOST_CHECK_EQUAL(result.evaluation, scalar_value_type::one());
        BOOST_REQUIRE_EQUAL(result.quotient.size(), 1);
        BOOST_CHECK(result.quotient[0].is_zero());
    }
}

BOOST_AUTO_TEST_SUITE_END()
