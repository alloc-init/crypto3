//---------------------------------------------------------------------------//
// Copyright (c) 2018-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
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

#define BOOST_TEST_MODULE sap_test

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>

#include <nil/crypto3/zk/snark/reductions/r1cs_to_sap.hpp>
#include <nil/crypto3/zk/snark/arithmetization/constraint_satisfaction_problems/r1cs.hpp>
#include <nil/crypto3/zk/snark/reductions/detail/sap_quotient.hpp>
#include <nil/crypto3/math/polynomial/backends/schoolbook_backend.hpp>
#include <nil/crypto3/math/polynomial/operations/lagrange_interpolation.hpp>

#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/algebra/curves/mnt4.hpp>
#include <nil/crypto3/algebra/fields/mnt4/base_field.hpp>
#include <nil/crypto3/algebra/fields/mnt4/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/mnt4.hpp>
#include <nil/crypto3/algebra/curves/params/multiexp/mnt4.hpp>
#include <nil/crypto3/algebra/curves/params/wnaf/mnt4.hpp>
#include <nil/crypto3/algebra/curves/mnt6.hpp>
#include <nil/crypto3/algebra/fields/mnt6/base_field.hpp>
#include <nil/crypto3/algebra/fields/mnt6/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/mnt6.hpp>
#include <nil/crypto3/algebra/curves/params/multiexp/mnt6.hpp>
#include <nil/crypto3/algebra/curves/params/wnaf/mnt6.hpp>

#include "../../systems/ppzksnark/r1cs_examples.hpp"

using namespace nil::crypto3::zk::snark;
using namespace nil::crypto3::algebra;

template<typename FieldType>
void test_sap(const std::size_t sap_degree, const std::size_t num_inputs, const bool binary_input) {
    /*
      We construct an instance where the SAP degree is <= sap_degree.
      The R1CS-to-SAP reduction produces SAPs with degree
        (2 * num_constraints + 2 * num_inputs + 1).
      So we generate an instance of R1CS where the number of constraints is
        (sap_degree - 1) / 2 - num_inputs.
    */
    // const std::size_t num_constraints = (sap_degree - 1) / 2 - num_inputs;
    const std::size_t num_constraints = 100;    // bc: reduced number of constraints so we can run this test regularly
    BOOST_CHECK(num_constraints >= 1);

    r1cs_example<FieldType> example;
    if (binary_input) {
        example = generate_r1cs_example_with_binary_input<FieldType>(num_constraints, num_inputs);
    } else {
        example = generate_r1cs_example_with_field_input<FieldType>(num_constraints, num_inputs);
    }
    BOOST_CHECK(example.constraint_system.is_satisfied(example.primary_input, example.auxiliary_input));

    const typename FieldType::value_type t = random_element<FieldType>(), d1 = random_element<FieldType>(),
                                         d2 = random_element<FieldType>();

    sap_instance<FieldType> sap_inst_1 = reductions::r1cs_to_sap<FieldType>::instance_map(example.constraint_system);

    sap_instance_evaluation<FieldType> sap_inst_2 =
        reductions::r1cs_to_sap<FieldType>::instance_map_with_evaluation(example.constraint_system, t);

    sap_witness<FieldType> sap_wit = reductions::r1cs_to_sap<FieldType>::witness_map(
        example.constraint_system, example.primary_input, example.auxiliary_input, d1, d2);

    BOOST_CHECK(sap_inst_1.is_satisfied(sap_wit));
    BOOST_CHECK(sap_inst_2.is_satisfied(sap_wit));
}

BOOST_AUTO_TEST_SUITE(sap_test_suite)

BOOST_AUTO_TEST_CASE(quotient_matches_polynomial_division) {
    using field_type = curves::mnt6<298>::scalar_field_type;
    using value_type = field_type::value_type;
    using polynomial_type = nil::crypto3::math::polynomial<value_type>;
    const nil::crypto3::math::polynomial_arithmetic::schoolbook_backend<value_type> reference;

    // Size 6 uses a step domain, where Z is not X^m - 1.
    for (const std::size_t m : {2, 4, 6, 8}) {
        const auto domain = nil::crypto3::math::make_evaluation_domain<field_type>(m);
        BOOST_REQUIRE(domain);
        BOOST_REQUIRE_EQUAL(domain->size(), m);
        const auto z = domain->get_vanishing_polynomial();
        for (const std::size_t nonzero_coefficients : {std::size_t(0), std::size_t(1), m}) {
            BOOST_TEST_CONTEXT("domain size " << m << ", nonzero coefficients " << nonzero_coefficients) {
                std::vector<value_type> a(m, value_type::zero());
                for (std::size_t i = 0; i < nonzero_coefficients; ++i) {
                    a[i] = value_type(i + 2);
                }
                polynomial_type a_squared, expected_h, c;
                reference.square(a_squared, polynomial_type(a));
                // Choosing C as the remainder makes A^2 - C exactly divisible by Z.
                nil::crypto3::math::division(expected_h, c, a_squared, z);
                c.resize(m, value_type::zero());
                expected_h.resize(m, value_type::zero());

                reductions::detail::compute_sap_quotient(*domain, a, c.get_storage());

                BOOST_REQUIRE_EQUAL(a.size(), m);
                BOOST_CHECK_EQUAL_COLLECTIONS(a.begin(), a.end(), expected_h.begin(), expected_h.end());
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(witness_quotient_matches_polynomial_division_with_blinding) {
    using field_type = curves::mnt6<298>::scalar_field_type;
    using value_type = field_type::value_type;
    using variable_type = nil::crypto3::math::linear_variable<field_type>;
    using polynomial_type = nil::crypto3::math::polynomial<value_type>;
    using reduction_type = reductions::r1cs_to_sap<field_type>;
    const nil::crypto3::math::polynomial_arithmetic::schoolbook_backend<value_type> reference;

    // x*x = y and y*x = z, with public x = 3 and private y = 9, z = 27.
    r1cs_constraint_system<field_type> cs;
    cs.primary_input_size = 1;
    cs.auxiliary_input_size = 2;
    cs.add_constraint({variable_type(1), variable_type(1), variable_type(2)});
    cs.add_constraint({variable_type(2), variable_type(1), variable_type(3)});
    const r1cs_primary_input<field_type> primary = {value_type(3)};
    const r1cs_auxiliary_input<field_type> auxiliary = {value_type(9), value_type(27)};
    BOOST_REQUIRE(cs.is_satisfied(primary, auxiliary));

    const auto domain = reduction_type::get_domain(cs);
    BOOST_REQUIRE_EQUAL(domain->size(), 8);
    // Four multiplication rows, three public-input rows, and one zero padding row.
    const std::array<int, 8> a_values = {6, 0, 12, 6, 1, 4, 2, 0};
    const std::array<int, 8> c_values = {36, 0, 144, 36, 1, 16, 4, 0};
    std::vector<std::pair<value_type, value_type>> a_points, c_points;
    for (std::size_t i = 0; i < domain->size(); ++i) {
        const auto x = domain->get_domain_element(i);
        a_points.emplace_back(x, value_type(a_values[i]));
        c_points.emplace_back(x, value_type(c_values[i]));
    }
    const auto unblinded_a = nil::crypto3::math::lagrange_interpolation(a_points);
    const auto unblinded_c = nil::crypto3::math::lagrange_interpolation(c_points);
    const auto z = domain->get_vanishing_polynomial();
    const std::vector<value_type> expected_assignment = {value_type(3), value_type(9),  value_type(27),
                                                         value_type(0), value_type(36), value_type(4)};
    const auto instance = reduction_type::instance_map_with_evaluation(cs, value_type(17));

    for (const auto &d1 : {value_type::zero(), value_type(5)}) {
        for (const auto &d2 : {value_type::zero(), value_type(7)}) {
            BOOST_TEST_CONTEXT("d1 " << d1 << ", d2 " << d2) {
                auto a = unblinded_a;
                auto c = unblinded_c;
                a.resize(domain->size() + 1, value_type::zero());
                c.resize(domain->size() + 1, value_type::zero());
                for (std::size_t i = 0; i < z.size(); ++i) {
                    a[i] += d1 * z[i];
                    c[i] += d2 * z[i];
                }
                polynomial_type numerator, expected_h, remainder;
                reference.square(numerator, a);
                numerator -= c;
                nil::crypto3::math::division(expected_h, remainder, numerator, z);
                BOOST_REQUIRE(remainder.is_zero());
                expected_h.resize(domain->size() + 1, value_type::zero());

                const auto witness = reduction_type::witness_map(cs, primary, auxiliary, d1, d2);

                BOOST_CHECK_EQUAL_COLLECTIONS(witness.coefficients_for_H.begin(), witness.coefficients_for_H.end(),
                                              expected_h.begin(), expected_h.end());
                BOOST_CHECK_EQUAL_COLLECTIONS(witness.coefficients_for_ACs.begin(), witness.coefficients_for_ACs.end(),
                                              expected_assignment.begin(), expected_assignment.end());
                BOOST_CHECK(instance.is_satisfied(witness));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(sap_test) {
    const std::size_t num_inputs = 10;

    /**
     * due to the specifics of our reduction, we can only get SAPs with odd
     * degrees, so we can only test "special" versions of the domains
     */

    using basic_curve_type = curves::mnt6<298>;

    const std::size_t basic_domain_size_special =
        (1ul << fields::arithmetic_params<basic_curve_type::scalar_field_type>::s) - 1ul;
    const std::size_t step_domain_size_special = (1ul << 10) + (1ul << 8) - 1ul;
    const std::size_t extended_domain_size_special =
        (1ul << (fields::arithmetic_params<basic_curve_type::scalar_field_type>::s + 1)) - 1ul;

    test_sap<typename basic_curve_type::scalar_field_type>(basic_domain_size_special, num_inputs, true);
    test_sap<typename basic_curve_type::scalar_field_type>(step_domain_size_special, num_inputs, true);
    test_sap<typename basic_curve_type::scalar_field_type>(extended_domain_size_special, num_inputs, true);

    test_sap<typename basic_curve_type::scalar_field_type>(basic_domain_size_special, num_inputs, false);
    test_sap<typename basic_curve_type::scalar_field_type>(step_domain_size_special, num_inputs, false);
    test_sap<typename basic_curve_type::scalar_field_type>(extended_domain_size_special, num_inputs, false);
}

BOOST_AUTO_TEST_SUITE_END()
