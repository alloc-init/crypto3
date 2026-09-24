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

#define BOOST_TEST_MODULE modified_sap_prover_test

#include <boost/multiprecision/integer.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/prover.hpp>

#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/operations/lagrange_interpolation.hpp>
#include <nil/crypto3/random/chacha_urbg.hpp>
#include <nil/crypto3/zk/algorithms/generate.hpp>
#include <nil/crypto3/zk/algorithms/prove.hpp>
#include <nil/crypto3/zk/algorithms/verify.hpp>
#include <nil/crypto3/zk/snark/arithmetization/bit_comparison.hpp>
#include <nil/crypto3/zk/snark/reductions/r1cs_to_modified_sap.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/generator.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/policy.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap/transcript.hpp>
#include <nil/crypto3/zk/snark/systems/ppsnark/modified_sap_snark.hpp>

namespace {

    using curve_type = nil::crypto3::algebra::curves::alt_bn128_254;
    using scalar_value_type = curve_type::scalar_field_type::value_type;
    using base_value_type = curve_type::base_field_type::value_type;
    using g1_value_type = curve_type::g1_type<>::value_type;
    using pairing_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_exact_pairing_policy;
    using transcript_policy_type = nil::crypto3::zk::snark::modified_sap_bn254_poseidon_transcript_policy;
    using policy_type =
        nil::crypto3::zk::snark::modified_sap_policy<curve_type, pairing_policy_type, transcript_policy_type>;
    using prover_type = nil::crypto3::zk::snark::modified_sap_prover<policy_type>;
    using scheme_type = nil::crypto3::zk::snark::modified_sap_snark<policy_type>;
    using deterministic_generator_type =
        nil::crypto3::zk::snark::detail::modified_sap_deterministic_generator<policy_type>;
    using reduction_type =
        nil::crypto3::zk::snark::reductions::modified_sap_to_polynomials<curve_type::scalar_field_type>;
    using r1cs_reduction_type =
        nil::crypto3::zk::snark::reductions::r1cs_to_modified_sap<curve_type::scalar_field_type>;
    using r1cs_system_type = nil::crypto3::zk::snark::r1cs_constraint_system<curve_type::scalar_field_type>;
    using r1cs_constraint_type = nil::crypto3::zk::snark::r1cs_constraint<curve_type::scalar_field_type>;
    using polynomial_type = reduction_type::polynomial_type;
    using auxiliary_input_type = policy_type::auxiliary_input_type;
    using proving_key_type = policy_type::proving_key_type;
    using proof_type = policy_type::proof_type;
    using system_type = policy_type::constraint_system_type;
    using constraint_type = system_type::constraint_type;
    using variable_type = constraint_type::variable_type;
    using combination_type = constraint_type::linear_combination_type;

    struct controlled_transcript_policy : transcript_policy_type {
        inline static challenge_type forced_challenge = challenge_type::zero();
        inline static std::size_t calls = 0;
        inline static verification_key_type last_verification_key;
        inline static commitment_type last_commitment = commitment_type::zero();
        inline static challenge_type last_public_input = challenge_type::zero();

        static challenge_type proof_challenge(const verification_key_type &verification_key,
                                              const commitment_type &P,
                                              const challenge_type &u) {
            ++calls;
            last_verification_key = verification_key;
            last_commitment = P;
            last_public_input = u;
            return forced_challenge;
        }
    };

    using controlled_policy_type =
        nil::crypto3::zk::snark::modified_sap_policy<curve_type, pairing_policy_type, controlled_transcript_policy>;
    using controlled_prover_type = nil::crypto3::zk::snark::modified_sap_prover<controlled_policy_type>;

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

    r1cs_system_type r1cs_test_system() {
        r1cs_system_type source;
        source.primary_input_size = 1;
        source.auxiliary_input_size = 2;
        // Source indices: 0 = one, 1 = u, 2 = x, 3 = y. Enforce (x + 1)*y = u.
        r1cs_constraint_type row {variable_type(0), variable_type(3), variable_type(1)};
        row.a.add_term(variable_type(2));
        source.add_constraint(row);
        return source;
    }

    proving_key_type generate_proving_key(const system_type &source) {
        // Fixed seed for reproducible tests with the production setup and transcript.
        const std::array<std::uint8_t, 32> seed = {};
        nil::crypto3::random::chacha_urbg<> random_source(seed);
        return scheme_type::generate(source, random_source).first;
    }

    deterministic_generator_type::trapdoor_type test_trapdoor() {
        return {scalar_value_type(2),
                scalar_value_type(3),
                {scalar_value_type(5), scalar_value_type(7), scalar_value_type(11), scalar_value_type(13)}};
    }

    scalar_value_type reference_quotient_at_tau(const polynomial_type &polynomial,
                                                const scalar_value_type &z,
                                                const scalar_value_type &tau) {
        // (X^j - z^j) / (X - z) = sum_(k=0..j-1) X^(j-1-k) * z^k, also valid at X = z.
        auto result = scalar_value_type::zero();
        for (std::size_t j = 1; j < polynomial.size(); ++j) {
            auto monomial_quotient = scalar_value_type::zero();
            for (std::size_t k = 0; k < j; ++k) {
                monomial_quotient += tau.pow(j - 1 - k) * z.pow(k);
            }
            result += polynomial[j] * monomial_quotient;
        }
        return result;
    }

}    // namespace

BOOST_AUTO_TEST_SUITE(modified_sap_prover_test_suite)

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

BOOST_AUTO_TEST_CASE(opening_commitment_matches_poseidon_and_scalar_reference) {
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const auto trapdoor = test_trapdoor();

    for (const std::size_t logical_rows : {1, 2, 3, 5}) {
        BOOST_TEST_CONTEXT("logical rows: " << logical_rows) {
            auto source = test_system();
            source.constraints.resize(logical_rows, binding_row());
            const auto key =
                deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor)
                    .first;
            const auto original_key = key;
            const auto committed = prover_type::commit(key, u, witness);
            const auto original_committed = committed;
            const auto z = transcript_policy_type::proof_challenge(key.verification_key, committed.P, u);
            const auto result = prover_type::open(key, u, committed);

            const std::array polynomials = {&committed.polynomials.A, &committed.polynomials.C,
                                            &committed.polynomials.H, &committed.Z};
            const std::array original_polynomials = {&original_committed.polynomials.A,
                                                     &original_committed.polynomials.C,
                                                     &original_committed.polynomials.H, &original_committed.Z};
            auto expected_scalar = scalar_value_type::zero();
            for (std::size_t i = 0; i < polynomials.size(); ++i) {
                BOOST_CHECK_EQUAL(result.evaluations[i], polynomials[i]->evaluate(z));
                expected_scalar += trapdoor.alpha[i] * reference_quotient_at_tau(*polynomials[i], z, trapdoor.tau);
                BOOST_CHECK(*polynomials[i] == *original_polynomials[i]);
            }
            // Q contains the alpha weights encoded in the T queries, with no gamma factor.
            BOOST_CHECK_EQUAL(result.Q, expected_scalar * g1_value_type::one());
            BOOST_CHECK(key == original_key);
            BOOST_CHECK_EQUAL(committed.P, original_committed.P);
        }
    }
}

BOOST_AUTO_TEST_CASE(opening_commitment_uses_transcript_inputs_and_accepts_edge_challenges) {
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const auto trapdoor = test_trapdoor();

    for (const std::size_t logical_rows : {2, 3, 5}) {
        auto source = test_system();
        source.constraints.resize(logical_rows, binding_row());
        const auto key =
            deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor)
                .first;
        const auto committed = controlled_prover_type::commit(key, u, witness);
        const auto domain = reduction_type::get_domain(source);
        std::vector<scalar_value_type> challenges {scalar_value_type::zero(), trapdoor.tau};
        for (std::size_t j = 0; j < domain->size(); ++j) {
            challenges.push_back(domain->get_domain_element(j));
        }
        const std::array polynomials = {&committed.polynomials.A, &committed.polynomials.C, &committed.polynomials.H,
                                        &committed.Z};

        for (const auto &z : challenges) {
            BOOST_TEST_CONTEXT("logical rows: " << logical_rows << ", challenge: " << z) {
                controlled_transcript_policy::forced_challenge = z;
                const auto calls_before = controlled_transcript_policy::calls;
                const auto result = controlled_prover_type::open(key, u, committed);
                BOOST_CHECK_EQUAL(controlled_transcript_policy::calls, calls_before + 1);
                BOOST_CHECK(controlled_transcript_policy::last_verification_key == key.verification_key);
                BOOST_CHECK_EQUAL(controlled_transcript_policy::last_commitment, committed.P);
                BOOST_CHECK_EQUAL(controlled_transcript_policy::last_public_input, u);

                auto expected_scalar = scalar_value_type::zero();
                for (std::size_t i = 0; i < polynomials.size(); ++i) {
                    BOOST_CHECK_EQUAL(result.evaluations[i], polynomials[i]->evaluate(z));
                    expected_scalar += trapdoor.alpha[i] * reference_quotient_at_tau(*polynomials[i], z, trapdoor.tau);
                }
                BOOST_CHECK_EQUAL(result.Q, expected_scalar * g1_value_type::one());
                if (domain->size() == 2) {
                    // H = 1: its canonical zero quotient must be skipped before accessing the empty T_H.
                    BOOST_REQUIRE(key.T_H.empty());
                    BOOST_CHECK_EQUAL(result.evaluations[2], scalar_value_type::one());
                    BOOST_CHECK_EQUAL(result.Q,
                                      (scalar_value_type(7) + scalar_value_type(13) * z) * g1_value_type::one());
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(opening_commitment_uses_query_prefixes_for_lower_degree_polynomials) {
    const auto u = scalar_value_type::one();
    const auto trapdoor = test_trapdoor();
    controlled_transcript_policy::forced_challenge = scalar_value_type::one();

    for (const std::size_t m : {4, 8}) {
        BOOST_TEST_CONTEXT("domain size: " << m) {
            system_type source {1, std::vector<constraint_type>(m, binding_row())};
            const auto domain = reduction_type::get_domain(source);
            // Prescribe A = X - 1 and C = A^2 - 1 = X^2 - 2X, giving H = 0.
            for (std::size_t j = 1; j < m; ++j) {
                const auto a = domain->get_domain_element(j) - scalar_value_type::one();
                source.constraints[j] = {combination_type(a), combination_type(a.squared() - u)};
            }
            const auto key =
                deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor)
                    .first;
            const auto committed = controlled_prover_type::commit(key, u, {u});
            BOOST_REQUIRE(committed.polynomials.H.is_zero());

            const auto result = controlled_prover_type::open(key, u, committed);
            // At z = 1, Q_A = 1, Q_C = X - 1, Q_H = 0 and Q_Z = 1 + X + ... + X^(m-1).
            const auto expected_scalar = trapdoor.alpha[0] + trapdoor.alpha[1] * (trapdoor.tau - u) +
                                         trapdoor.alpha[3] * (trapdoor.tau.pow(m) - u) * (trapdoor.tau - u).inversed();
            BOOST_CHECK_EQUAL(result.Q, expected_scalar * g1_value_type::one());
            BOOST_CHECK(result.evaluations[0].is_zero());
            BOOST_CHECK_EQUAL(result.evaluations[1], -u);
            BOOST_CHECK(result.evaluations[2].is_zero());
            BOOST_CHECK(result.evaluations[3].is_zero());
        }
    }
}

BOOST_AUTO_TEST_CASE(opening_commitment_accepts_identity_bases_and_cancellation) {
    auto source = test_system();
    source.constraints.resize(2);
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    controlled_transcript_policy::forced_challenge = scalar_value_type::zero();

    for (const bool cancellation : {false, true}) {
        BOOST_TEST_CONTEXT("cancellation: " << cancellation) {
            auto trapdoor = test_trapdoor();
            if (cancellation) {
                // At z = 0, Q_A = -1, Q_C = -2 and Q_Z(tau) = tau, so all three contributions cancel.
                trapdoor.alpha[3] =
                    (trapdoor.alpha[0] + scalar_value_type(2) * trapdoor.alpha[1]) * trapdoor.tau.inversed();
            } else {
                trapdoor.alpha.fill(scalar_value_type::zero());
            }
            const auto key =
                deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor)
                    .first;
            for (const auto *query : {&key.T_A, &key.T_C, &key.T_Z}) {
                for (const auto &point : *query) {
                    BOOST_REQUIRE(point.is_zero() == !cancellation);
                }
            }
            const auto committed = controlled_prover_type::commit(key, u, witness);
            const auto result = controlled_prover_type::open(key, u, committed);
            BOOST_CHECK(result.Q.is_zero());
        }
    }
}

BOOST_AUTO_TEST_CASE(opening_commitment_rejects_quotients_exceeding_query_lengths) {
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    controlled_transcript_policy::forced_challenge = scalar_value_type::zero();
    const std::array<const char *, 4> names = {"A", "C", "H", "Z"};

    for (const std::size_t logical_rows : {2, 3}) {
        auto source = test_system();
        source.constraints.resize(logical_rows);
        const auto key = generate_proving_key(source);
        const auto committed = controlled_prover_type::commit(key, u, witness);
        const std::array queries = {&key.T_A, &key.T_C, &key.T_H, &key.T_Z};
        for (std::size_t i = 0; i < queries.size(); ++i) {
            BOOST_TEST_CONTEXT("logical rows: " << logical_rows << ", polynomial: " << names[i]) {
                auto changed = committed;
                const std::array polynomials = {&changed.polynomials.A, &changed.polynomials.C, &changed.polynomials.H,
                                                &changed.Z};
                // At z = 0 this monomial has a nonzero quotient one coefficient longer than its query.
                *polynomials[i] = polynomial_type(queries[i]->size() + 2, scalar_value_type::zero());
                polynomials[i]->back() = scalar_value_type::one();
                BOOST_CHECK_THROW(controlled_prover_type::open(key, u, changed), std::invalid_argument);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(opening_commitment_rejects_shortened_queries) {
    const auto original = generate_proving_key(test_system());
    const scalar_value_type u(3);
    const auto committed = controlled_prover_type::commit(original, u, {u, scalar_value_type(2), scalar_value_type(1)});
    controlled_transcript_policy::forced_challenge = scalar_value_type::zero();
    const std::array queries = {std::pair {"T_A", &proving_key_type::T_A}, std::pair {"T_C", &proving_key_type::T_C},
                                std::pair {"T_H", &proving_key_type::T_H}, std::pair {"T_Z", &proving_key_type::T_Z}};
    for (const auto &[name, member] : queries) {
        BOOST_TEST_CONTEXT("query: " << name) {
            auto key = original;
            (key.*member).clear();
            BOOST_CHECK_THROW(controlled_prover_type::open(key, u, committed), std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(process_matches_scalar_reference) {
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const auto trapdoor = test_trapdoor();

    for (const std::size_t logical_rows : {1, 2, 3, 5}) {
        BOOST_TEST_CONTEXT("logical rows: " << logical_rows) {
            auto source = test_system();
            source.constraints.resize(logical_rows, binding_row());
            const auto key =
                deterministic_generator_type::process(source, transcript_policy_type::circuit_digest(source), trapdoor)
                    .first;
            const auto proof = prover_type::process(key, u, witness);

            // Derive expected commitments from scalar evaluations, without commit(), open() or MSMs.
            const auto polynomials = reduction_type::witness_map(source, u, witness);
            const auto Z = reduction_type::get_domain(source)->get_vanishing_polynomial();
            const auto p_scalar = trapdoor.gamma * (trapdoor.alpha[0] * polynomials.A.evaluate(trapdoor.tau) +
                                                    trapdoor.alpha[1] * polynomials.C.evaluate(trapdoor.tau) +
                                                    trapdoor.alpha[2] * polynomials.H.evaluate(trapdoor.tau));
            const auto expected_P = p_scalar * g1_value_type::one();
            const auto z = transcript_policy_type::proof_challenge(key.verification_key, expected_P, u);
            auto q_scalar = scalar_value_type::zero();
            const std::array opened_polynomials = {&polynomials.A, &polynomials.C, &polynomials.H, &Z};
            for (std::size_t i = 0; i < opened_polynomials.size(); ++i) {
                q_scalar += trapdoor.alpha[i] * reference_quotient_at_tau(*opened_polynomials[i], z, trapdoor.tau);
            }

            BOOST_CHECK_EQUAL(proof.P, expected_P);
            BOOST_CHECK_EQUAL(proof.Q, q_scalar * g1_value_type::one());
            BOOST_CHECK_EQUAL(proof.v_A, polynomials.A.evaluate(z));
            BOOST_CHECK_EQUAL(proof.v_C, polynomials.C.evaluate(z));
            BOOST_CHECK_EQUAL(proof.v_H, polynomials.H.evaluate(z));
            BOOST_CHECK_EQUAL(proof.v_Z, Z.evaluate(z));
        }
    }
}

BOOST_AUTO_TEST_CASE(prove_entry_points_are_deterministic_and_preserve_inputs) {
    const auto key = generate_proving_key(test_system());
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const auto original_key = key;
    const auto original_u = u;
    const auto original_witness = witness;

    const auto proof = prover_type::process(key, u, witness);
    BOOST_CHECK(prover_type::process(key, u, witness) == proof);
    BOOST_CHECK(scheme_type::prove(key, u, witness) == proof);
    BOOST_CHECK(nil::crypto3::zk::prove<scheme_type>(key, u, witness) == proof);
    BOOST_CHECK(key == original_key);
    BOOST_CHECK_EQUAL(u, original_u);
    BOOST_CHECK(witness == original_witness);
}

BOOST_AUTO_TEST_CASE(prove_entry_points_reject_invalid_inputs_without_replacing_a_proof) {
    const auto key = generate_proving_key(test_system());
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const auto original_proof = prover_type::process(key, u, witness);

    const auto expect_rejection = [&](const proving_key_type &supplied_key,
                                      const scalar_value_type &supplied_u,
                                      const auxiliary_input_type &supplied_witness) {
        const auto original_key = supplied_key;
        const auto original_u = supplied_u;
        const auto original_witness = supplied_witness;
        auto proof = original_proof;
        BOOST_CHECK_THROW(proof = prover_type::process(supplied_key, supplied_u, supplied_witness),
                          std::invalid_argument);
        BOOST_CHECK(proof == original_proof);
        BOOST_CHECK_THROW(proof = scheme_type::prove(supplied_key, supplied_u, supplied_witness),
                          std::invalid_argument);
        BOOST_CHECK(proof == original_proof);
        BOOST_CHECK_THROW(proof = nil::crypto3::zk::prove<scheme_type>(supplied_key, supplied_u, supplied_witness),
                          std::invalid_argument);
        BOOST_CHECK(proof == original_proof);
        BOOST_CHECK(supplied_key == original_key);
        BOOST_CHECK_EQUAL(supplied_u, original_u);
        BOOST_CHECK(supplied_witness == original_witness);
    };

    auto shortened_key = key;
    shortened_key.T_Z.pop_back();
    expect_rejection(shortened_key, u, witness);

    auto changed_digest = key;
    changed_digest.verification_key.circuit_digest += base_value_type::one();
    expect_rejection(changed_digest, u, witness);

    auto short_witness = witness;
    short_witness.pop_back();
    expect_rejection(key, u, short_witness);

    auto wrong_binding = witness;
    wrong_binding[0] = scalar_value_type(5);
    expect_rejection(key, u, wrong_binding);

    auto unsatisfied = witness;
    unsatisfied[1] += scalar_value_type::one();
    expect_rejection(key, u, unsatisfied);

    const auto binding_key = generate_proving_key({1, {binding_row()}});
    const scalar_value_type even_u(2);
    expect_rejection(binding_key, even_u, {even_u});
}

BOOST_AUTO_TEST_CASE(prove_propagates_opening_failure_without_replacing_a_proof) {
    struct throwing_transcript_policy : transcript_policy_type {
        static challenge_type proof_challenge(const verification_key_type &, const commitment_type &,
                                              const challenge_type &) {
            throw std::runtime_error("test transcript failure");
        }
    };
    using throwing_policy_type =
        nil::crypto3::zk::snark::modified_sap_policy<curve_type, pairing_policy_type, throwing_transcript_policy>;
    using throwing_scheme_type = nil::crypto3::zk::snark::modified_sap_snark<throwing_policy_type>;

    const auto key = generate_proving_key(test_system());
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const auto original_key = key;
    const auto original_u = u;
    const auto original_witness = witness;
    const auto original_proof = prover_type::process(key, u, witness);
    auto proof = original_proof;

    BOOST_CHECK_THROW(proof = nil::crypto3::zk::prove<throwing_scheme_type>(key, u, witness), std::runtime_error);
    BOOST_CHECK(proof == original_proof);
    BOOST_CHECK(key == original_key);
    BOOST_CHECK_EQUAL(u, original_u);
    BOOST_CHECK(witness == original_witness);
}

BOOST_AUTO_TEST_CASE(generated_proofs_verify_across_domain_sizes_and_setup_seeds) {
    const scalar_value_type u(3);
    const std::array sizes = {std::pair {1, 2}, std::pair {2, 2}, std::pair {3, 4}, std::pair {4, 4},
                              std::pair {5, 8}, std::pair {8, 8}, std::pair {9, 16}};
    for (const unsigned seed_tag : {0, 1}) {
        // Reproducible streams exercise the production randomized setup.
        std::array<std::uint8_t, 32> seed = {};
        seed[0] = static_cast<std::uint8_t>(seed_tag);
        nil::crypto3::random::chacha_urbg<> random_source(seed);
        for (const auto &[logical_rows, domain_size] : sizes) {
            BOOST_TEST_CONTEXT("seed: " << seed_tag << ", logical rows: " << logical_rows) {
                auto source = test_system();
                source.constraints.resize(logical_rows, binding_row());
                auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
                if (logical_rows == 1) {
                    // The minimum circuit has only the public-input binding row and witness entry.
                    source.witness_size = 1;
                    witness.resize(1);
                }
                const auto keys = nil::crypto3::zk::generate<scheme_type>(source, random_source);
                BOOST_REQUIRE_EQUAL(keys.second.domain_size, domain_size);
                if (domain_size == 2) {
                    BOOST_REQUIRE(keys.first.T_H.empty());
                }

                const auto proof = nil::crypto3::zk::prove<scheme_type>(keys.first, u, witness);
                BOOST_REQUIRE(scheme_type::verify(keys.second, u, proof));
                if (seed_tag == 0 && logical_rows == 1) {
                    // Exercise the generic verification adapter once.
                    BOOST_CHECK(nil::crypto3::zk::verify<scheme_type>(keys.second, u, proof));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(generated_proofs_verify_with_zero_private_witness_and_zero_h) {
    system_type source {2, {binding_row(), {combination_type(variable_type(1)), combination_type(-variable_type(0))}}};
    source.constraints.resize(5, binding_row());
    const auto key = generate_proving_key(source);

    for (const auto &u : {scalar_value_type(1), scalar_value_type(3), scalar_value_type(5), -scalar_value_type(2)}) {
        BOOST_TEST_CONTEXT("public input: " << u) {
            // The canonical representative of -2 is r - 2, which is odd.
            const auto proof = scheme_type::prove(key, u, {u, scalar_value_type::zero()});
            // Every row has A = 0 and C = -u, so the quotient H is zero.
            BOOST_CHECK(proof.v_A.is_zero());
            BOOST_CHECK_EQUAL(proof.v_C, -u);
            BOOST_CHECK(proof.v_H.is_zero());
            BOOST_CHECK(scheme_type::verify(key.verification_key, u, proof));
        }
    }
}

BOOST_AUTO_TEST_CASE(generated_proofs_reject_changed_commitments_and_evaluations) {
    const auto key = generate_proving_key(test_system());
    const scalar_value_type u(3);
    const auto proof = scheme_type::prove(key, u, {u, scalar_value_type(2), scalar_value_type(1)});
    BOOST_REQUIRE(scheme_type::verify(key.verification_key, u, proof));

    const std::array commitments = {std::pair {"P", &proof_type::P}, std::pair {"Q", &proof_type::Q}};
    for (const auto &[name, member] : commitments) {
        BOOST_TEST_CONTEXT("commitment: " << name) {
            auto changed = proof;
            changed.*member += g1_value_type::one();
            BOOST_CHECK(!scheme_type::verify(key.verification_key, u, changed));
        }
    }
    const std::array evaluations = {std::pair {"v_A", &proof_type::v_A}, std::pair {"v_C", &proof_type::v_C},
                                    std::pair {"v_H", &proof_type::v_H}, std::pair {"v_Z", &proof_type::v_Z}};
    for (const auto &[name, member] : evaluations) {
        BOOST_TEST_CONTEXT("evaluation: " << name) {
            auto changed = proof;
            changed.*member += scalar_value_type::one();
            BOOST_CHECK(!scheme_type::verify(key.verification_key, u, changed));
        }
    }

    // Changing C and H together preserves the arithmetic equation but invalidates the opening commitment.
    auto changed = proof;
    changed.v_H += scalar_value_type::one();
    changed.v_C -= proof.v_Z;
    BOOST_REQUIRE_EQUAL(changed.v_A.squared() - changed.v_C, changed.v_H * changed.v_Z + u);
    BOOST_CHECK(!scheme_type::verify(key.verification_key, u, changed));
}

BOOST_AUTO_TEST_CASE(generated_proofs_are_bound_to_the_expected_public_input) {
    const auto key = generate_proving_key(test_system());
    const scalar_value_type u(3);
    const auto proof = scheme_type::prove(key, u, {u, scalar_value_type(2), scalar_value_type(1)});
    BOOST_REQUIRE(scheme_type::verify(key.verification_key, u, proof));

    for (const auto &changed_u : {scalar_value_type(1), scalar_value_type(5)}) {
        BOOST_TEST_CONTEXT("changed public input: " << changed_u) {
            BOOST_CHECK(!scheme_type::verify(key.verification_key, changed_u, proof));

            // Restore the arithmetic equation for the new public input; the transcript and openings still bind u.
            auto changed = proof;
            changed.v_C += u - changed_u;
            BOOST_REQUIRE_EQUAL(changed.v_A.squared() - changed.v_C, changed.v_H * changed.v_Z + changed_u);
            BOOST_CHECK(!scheme_type::verify(key.verification_key, changed_u, changed));
        }
    }
}

BOOST_AUTO_TEST_CASE(generated_proofs_are_bound_to_the_circuit_and_setup) {
    const auto source = test_system();
    const scalar_value_type u(3);
    const auxiliary_input_type witness = {u, scalar_value_type(2), scalar_value_type(1)};
    const std::array<std::uint8_t, 32> seed = {};
    nil::crypto3::random::chacha_urbg<> random_source(seed);
    const auto keys = scheme_type::generate(source, random_source);
    const auto proof = scheme_type::prove(keys.first, u, witness);
    BOOST_REQUIRE(scheme_type::verify(keys.second, u, proof));

    // A second setup for the same circuit uses fresh draws from the ChaCha stream.
    const auto fresh_keys = scheme_type::generate(source, random_source);
    const auto fresh_proof = scheme_type::prove(fresh_keys.first, u, witness);
    BOOST_REQUIRE(scheme_type::verify(fresh_keys.second, u, fresh_proof));
    BOOST_REQUIRE(!(fresh_keys.second == keys.second));
    BOOST_CHECK(!scheme_type::verify(fresh_keys.second, u, proof));
    BOOST_CHECK(!scheme_type::verify(keys.second, u, fresh_proof));

    // Negating A in one row preserves satisfaction but changes the circuit digest.
    auto changed_source = source;
    changed_source.constraints[1].a = combination_type(-variable_type(1));
    nil::crypto3::random::chacha_urbg<> repeated_random_source(seed);
    const auto changed_keys = scheme_type::generate(changed_source, repeated_random_source);
    const auto changed_proof = scheme_type::prove(changed_keys.first, u, witness);
    BOOST_REQUIRE(scheme_type::verify(changed_keys.second, u, changed_proof));
    BOOST_REQUIRE_NE(changed_keys.second.circuit_digest, keys.second.circuit_digest);

    // Equal dimensions and identical setup draws make the verification keys differ only in their circuit digest.
    auto matching_parameters = changed_keys.second;
    matching_parameters.circuit_digest = keys.second.circuit_digest;
    BOOST_REQUIRE(matching_parameters == keys.second);
    BOOST_CHECK(!scheme_type::verify(changed_keys.second, u, proof));
    BOOST_CHECK(!scheme_type::verify(keys.second, u, changed_proof));
}

BOOST_AUTO_TEST_CASE(r1cs_frontend_reuses_setup_for_multiple_assignments) {
    // One public input u and two private variables x, y, constrained by (x + 1)*y = u.
    const auto source = r1cs_test_system();
    const auto converted = r1cs_reduction_type::instance_map(source);
    const auto key = generate_proving_key(converted);
    BOOST_REQUIRE_EQUAL(key.constraint_system.num_constraints(), 564);
    BOOST_REQUIRE_EQUAL(key.verification_key.num_variables, 510);
    BOOST_REQUIRE_EQUAL(key.verification_key.domain_size, 1024);

    // Equivalent unsorted, repeated, zero and cancelling terms use the same converted circuit and key.
    auto raw_source = source;
    raw_source.constraints[0].a.add_term(variable_type(2), scalar_value_type(2));
    raw_source.constraints[0].a.add_term(variable_type(0), scalar_value_type::zero());
    raw_source.constraints[0].a.add_term(variable_type(2), -scalar_value_type(2));
    BOOST_REQUIRE(r1cs_reduction_type::instance_map(raw_source) == converted);

    for (const auto &u : {scalar_value_type(1), scalar_value_type(3), scalar_value_type(5), -scalar_value_type(2)}) {
        BOOST_TEST_CONTEXT("public input: " << u) {
            const auto witness = r1cs_reduction_type::witness_map(
                raw_source, {u}, {u - scalar_value_type::one(), scalar_value_type::one()});
            const auto proof = nil::crypto3::zk::prove<scheme_type>(key, u, witness);
            BOOST_REQUIRE(nil::crypto3::zk::verify<scheme_type>(key.verification_key, u, proof));
            const auto changed_u = u == scalar_value_type(3) ? scalar_value_type(5) : scalar_value_type(3);
            BOOST_CHECK(!scheme_type::verify(key.verification_key, changed_u, proof));
        }
    }
}

BOOST_AUTO_TEST_CASE(r1cs_frontend_prover_rejects_changed_witness_entries) {
    const auto source = r1cs_test_system();
    const auto key = generate_proving_key(r1cs_reduction_type::instance_map(source));
    const scalar_value_type u(3);
    const auto witness = r1cs_reduction_type::witness_map(source, {u}, {scalar_value_type(2), scalar_value_type(1)});
    const auto proof = scheme_type::prove(key, u, witness);
    BOOST_REQUIRE(scheme_type::verify(key.verification_key, u, proof));

    // N = 3; recovery has 254 bits, 99 new markers and 153 comparison auxiliaries.
    const std::array changes = {
        std::pair {"public-input entry", std::size_t(0)},     std::pair {"source variable", std::size_t(1)},
        std::pair {"recovered one", std::size_t(3)},          std::pair {"comparison marker", std::size_t(257)},
        std::pair {"comparison auxiliary", std::size_t(356)}, std::pair {"source auxiliary", std::size_t(509)}};
    BOOST_REQUIRE_EQUAL(witness.size(), 510);
    for (const auto &[name, index] : changes) {
        BOOST_TEST_CONTEXT("changed entry: " << name) {
            auto changed = witness;
            // Keep the recovered bit Boolean when changing the source's constant-one wire.
            changed[index] = index == 3 ? scalar_value_type::zero() : changed[index] + scalar_value_type::one();
            BOOST_CHECK_THROW(scheme_type::prove(key, u, changed), std::invalid_argument);
        }
    }
    BOOST_CHECK_THROW(scheme_type::prove(key, scalar_value_type(5), witness), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(r1cs_frontend_proofs_cover_empty_sources_and_domain_boundaries) {
    const scalar_value_type u(3);
    const std::array sizes = {std::pair {0, 1024}, std::pair {231, 1024}, std::pair {232, 2048}};
    for (const auto &[source_rows, domain_size] : sizes) {
        BOOST_TEST_CONTEXT("source rows: " << source_rows << ", domain size: " << domain_size) {
            auto source = r1cs_test_system();
            const auto multiplication = source.constraints.front();
            source.constraints.resize(source_rows, multiplication);
            const auto converted = r1cs_reduction_type::instance_map(source);
            BOOST_REQUIRE_EQUAL(converted.num_constraints(), 562 + 2 * source_rows);
            const auto key = generate_proving_key(converted);
            BOOST_REQUIRE_EQUAL(key.verification_key.domain_size, domain_size);
            // Setup keeps the logical rows; polynomial construction supplies binding-row padding.
            BOOST_CHECK_EQUAL(key.constraint_system.num_constraints(), converted.num_constraints());

            const auto witness =
                r1cs_reduction_type::witness_map(source, {u}, {scalar_value_type(2), scalar_value_type(1)});
            const auto proof = scheme_type::prove(key, u, witness);
            BOOST_CHECK(scheme_type::verify(key.verification_key, u, proof));
        }
    }
}

BOOST_AUTO_TEST_CASE(r1cs_bitcoin_amount_range_proofs) {
    // Bitcoin amounts satisfy 0 <= amount <= MAX_MONEY, measured in satoshis.
    constexpr std::uint64_t satoshis_per_bitcoin = 100'000'000;
    constexpr std::uint64_t max_money = 21'000'000 * satoshis_per_bitcoin;
    constexpr std::uint64_t bound = max_money + 1;
    const std::size_t bit_count = boost::multiprecision::msb(bound) + 1;
    BOOST_REQUIRE_EQUAL(bit_count, 51);

    r1cs_system_type source;
    source.primary_input_size = 1;
    source.auxiliary_input_size = 1 + bit_count;
    // Source indices: 0 = one, 1 = u, 2 = amount, then little-endian bits and comparison markers.
    const variable_type amount_variable(2);
    constexpr std::size_t first_bit = 3;
    for (std::size_t bit = 0; bit < bit_count; ++bit) {
        // b*(b - 1) = 0.
        r1cs_constraint_type boolean;
        boolean.a.add_term(variable_type(first_bit + bit));
        boolean.b.add_term(variable_type(0), -scalar_value_type::one());
        boolean.b.add_term(variable_type(first_bit + bit));
        source.add_constraint(boolean);
    }

    // The 51-bit sum is smaller than the scalar-field modulus, so packing cannot wrap.
    r1cs_constraint_type packing;
    packing.b.add_term(variable_type(0));
    packing.c.add_term(amount_variable);
    auto weight = scalar_value_type::one();
    for (std::size_t bit = 0; bit < bit_count; ++bit) {
        packing.a.add_term(variable_type(first_bit + bit), weight);
        weight += weight;
    }
    source.add_constraint(packing);

    // u = 2*amount + 1 satisfies the inner SNARK's odd-public-input rule for every valid amount.
    r1cs_constraint_type public_binding {variable_type(0), variable_type(0), variable_type(1)};
    public_binding.a.add_term(amount_variable, scalar_value_type(2));
    source.add_constraint(public_binding);

    // The shared comparator enforces amount < MAX_MONEY + 1, including the upper endpoint.
    variable_type marker(0);
    nil::crypto3::zk::snark::for_each_bit_comparison(
        bound, bit_count, curve_type::scalar_field_type::modulus,
        [&](std::size_t bit) { marker = variable_type(first_bit + bit); },
        [&](std::size_t bit) {
            const variable_type next_marker(source.num_variables() + 1);
            ++source.auxiliary_input_size;
            source.add_constraint({marker, variable_type(first_bit + bit), next_marker});
            marker = next_marker;
        },
        [&](std::size_t begin, std::size_t end) {
            r1cs_constraint_type zero_product;
            zero_product.a.add_term(marker);
            for (std::size_t bit = begin; bit < end; ++bit) {
                zero_product.b.add_term(variable_type(first_bit + bit));
            }
            source.add_constraint(zero_product);
        });
    BOOST_REQUIRE(source.is_valid());

    // 1 public input + 1 amount + 51 bits + 19 comparison markers; the implicit one is excluded.
    BOOST_REQUIRE_EQUAL(source.num_variables(), 72);
    BOOST_REQUIRE_EQUAL(source.auxiliary_input_size, 71);
    // 51 Booleanity + 1 packing + 1 public binding + 29 comparison constraints.
    BOOST_REQUIRE_EQUAL(source.num_constraints(), 82);

    const auto make_auxiliary = [&](const scalar_value_type &amount) {
        std::vector<scalar_value_type> auxiliary = {amount};
        auxiliary.reserve(source.auxiliary_input_size);
        const auto integer = amount.to_integral();
        for (std::size_t bit = 0; bit < bit_count; ++bit) {
            auxiliary.emplace_back(boost::multiprecision::bit_test(integer, bit));
        }
        auto marker_value = scalar_value_type::one();
        nil::crypto3::zk::snark::for_each_bit_comparison(
            bound, bit_count, curve_type::scalar_field_type::modulus,
            [&](std::size_t bit) { marker_value = auxiliary[1 + bit]; },
            [&](std::size_t bit) {
                marker_value *= auxiliary[1 + bit];
                auxiliary.push_back(marker_value);
            },
            // No range check here: the circuit must reject out-of-range assignments itself.
            [](std::size_t, std::size_t) { });
        return auxiliary;
    };

    // One public-API setup serves zero, one satoshi, one BTC and both upper-bound cases.
    const auto converted = r1cs_reduction_type::instance_map(source);
    // 72 source values + 506 constant-recovery values + 82 multiplication auxiliaries.
    BOOST_REQUIRE_EQUAL(converted.num_variables(), 660);
    // 1 public binding + 561 constant-recovery rows + 2*82 source rows, before padding.
    BOOST_REQUIRE_EQUAL(converted.num_constraints(), 726);
    const auto key = generate_proving_key(converted);
    BOOST_REQUIRE_EQUAL(key.verification_key.domain_size, 1024);
    const std::array<std::uint64_t, 5> valid_amounts = {0, 1, satoshis_per_bitcoin, max_money - 1, max_money};
    for (const auto amount : valid_amounts) {
        BOOST_TEST_CONTEXT("satoshis: " << amount) {
            const auto u = scalar_value_type(2 * amount + 1);
            const auto auxiliary = make_auxiliary(scalar_value_type(amount));
            BOOST_REQUIRE(source.is_satisfied({u}, auxiliary));
            const auto witness = r1cs_reduction_type::witness_map(source, {u}, auxiliary);
            BOOST_REQUIRE(converted.is_satisfied(u, witness));
            const auto proof = scheme_type::prove(key, u, witness);
            BOOST_REQUIRE(scheme_type::verify(key.verification_key, u, proof));
            BOOST_CHECK(!scheme_type::verify(key.verification_key, u + scalar_value_type(2), proof));
        }
    }

    // MAX_MONEY + 1 and 2^51 - 1 fit the bit width: rejection must come from the comparison.
    auto without_comparison = source;
    without_comparison.constraints.resize(bit_count + 2);
    const std::uint64_t above_bit_width = std::uint64_t(1) << bit_count;
    for (const auto amount : {max_money + 1, above_bit_width - 1}) {
        const auto u = scalar_value_type(2 * amount + 1);
        BOOST_REQUIRE(without_comparison.is_satisfied({u}, make_auxiliary(scalar_value_type(amount))));
    }
    const std::array invalid_amounts = {-scalar_value_type::one(), scalar_value_type(max_money + 1),
                                        scalar_value_type(above_bit_width - 1), scalar_value_type(above_bit_width)};
    for (const auto &amount : invalid_amounts) {
        BOOST_TEST_CONTEXT("out-of-range amount: " << amount) {
            const auto u = scalar_value_type(2) * amount + scalar_value_type::one();
            const auto auxiliary = make_auxiliary(amount);
            BOOST_CHECK(!source.is_satisfied({u}, auxiliary));
            BOOST_CHECK_THROW(r1cs_reduction_type::witness_map(source, {u}, auxiliary), std::invalid_argument);
        }
    }

    // Replace bits [0, 1] by [2, 0]: the packed amount stays 2, but Booleanity must fail.
    const scalar_value_type u(5);
    auto auxiliary = make_auxiliary(scalar_value_type(2));
    auto witness = r1cs_reduction_type::witness_map(source, {u}, auxiliary);
    auxiliary[1] = scalar_value_type(2);
    auxiliary[2] = scalar_value_type::zero();
    BOOST_CHECK(!source.is_satisfied({u}, auxiliary));
    BOOST_CHECK_THROW(r1cs_reduction_type::witness_map(source, {u}, auxiliary), std::invalid_argument);
    // Also reject a modified-SAP witness supplied directly to the public prover.
    witness[first_bit - 1] = auxiliary[1];
    witness[first_bit] = auxiliary[2];
    BOOST_CHECK_THROW(scheme_type::prove(key, u, witness), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
