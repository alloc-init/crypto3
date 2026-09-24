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

#define BOOST_TEST_MODULE r1cs_to_modified_sap_test

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/zk/snark/reductions/r1cs_to_modified_sap.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;
    namespace snark = nil::crypto3::zk::snark;
    using field_type = fields::alt_bn128_scalar_field<254>;
    using value_type = field_type::value_type;
    using recovery_type = snark::reductions::detail::modified_sap_constant_recovery<field_type>;

    // Use Crypto3's prime-field arithmetic with small moduli for exhaustive checks.
    template<unsigned Modulus>
    class test_prime_field : public fields::field<8> {
    public:
        using integral_type = fields::field<8>::integral_type;
        constexpr static integral_type modulus = Modulus;
        constexpr static integral_type group_order_minus_one_half = (modulus - 1u) / 2;
        constexpr static modular_params_type modulus_params = modulus.backend();
        using modular_type = boost::multiprecision::number<boost::multiprecision::backends::modular_adaptor<
            modular_backend, boost::multiprecision::backends::modular_params_ct<modular_backend, modulus_params>>>;
        using value_type = fields::detail::element_fp<fields::params<test_prime_field>>;
    };

    template<typename FieldType>
    snark::modified_sap_constraint_system<FieldType> recovery_system(std::size_t first_bit) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        using variable = typename recovery::variable_type;
        snark::modified_sap_constraint_system<FieldType> cs;
        cs.witness_size = first_bit + recovery::witness_size;
        cs.constraints.push_back({{}, -variable(0)});
        recovery::append_constraints(cs.constraints, first_bit);
        return cs;
    }

    template<typename FieldType>
    bool rows_satisfied(const snark::modified_sap_constraint_system<FieldType> &cs,
                        const typename FieldType::value_type &u,
                        const std::vector<typename FieldType::value_type> &witness, std::size_t begin,
                        std::size_t end) {
        for (std::size_t row = begin; row < end; ++row) {
            if (cs.constraints[row].a.evaluate(witness).squared() - cs.constraints[row].c.evaluate(witness) != u) {
                return false;
            }
        }
        return true;
    }

    template<typename FieldType>
    void complete_product_auxiliaries(const snark::modified_sap_constraint_system<FieldType> &cs,
                                      std::vector<typename FieldType::value_type> &witness,
                                      std::size_t first_bit) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        std::size_t auxiliary = first_bit + recovery::bit_count + recovery::marker_count;
        // Every second product row forces t = A(w)^2 when w[0] = u. Satisfy it even for forged bits/markers.
        for (std::size_t row = recovery::bit_count + 3; row < cs.constraints.size(); row += 2) {
            witness[auxiliary++] = cs.constraints[row].a.evaluate(witness).squared();
        }
    }

    template<typename FieldType>
    void complete_comparison_witness(const snark::modified_sap_constraint_system<FieldType> &cs,
                                     std::vector<typename FieldType::value_type> &witness,
                                     std::size_t first_bit) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        auto marker = witness[first_bit + recovery::bit_count - 1];
        std::size_t next_marker = first_bit + recovery::bit_count;
        // Reference markers are prefix products over the modulus's internal one bits.
        for (std::size_t bit = recovery::bit_count - 1; bit > 1;) {
            --bit;
            if (boost::multiprecision::bit_test(FieldType::modulus, bit)) {
                marker *= witness[first_bit + bit];
                witness[next_marker++] = marker;
            }
        }
        complete_product_auxiliaries(cs, witness, first_bit);
    }

    template<typename FieldType>
    std::vector<typename FieldType::value_type>
        witness_from_integer(const snark::modified_sap_constraint_system<FieldType> &cs,
                             const typename FieldType::integral_type &integer, std::size_t first_bit) {
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<FieldType>;
        using value = typename FieldType::value_type;
        std::vector<value> witness(cs.witness_size, value::zero());
        witness[0] = value(integer);
        for (std::size_t bit = 0; bit < recovery::bit_count; ++bit) {
            witness[first_bit + bit] = value(boost::multiprecision::bit_test(integer, bit));
        }
        complete_comparison_witness(cs, witness, first_bit);
        return witness;
    }

    template<unsigned Modulus>
    void test_all_bit_patterns() {
        using field = test_prime_field<Modulus>;
        using value = typename field::value_type;
        using recovery = snark::reductions::detail::modified_sap_constant_recovery<field>;
        const auto cs = recovery_system<field>(1);
        BOOST_REQUIRE(cs.is_valid());
        for (unsigned integer = 0; integer < (1u << recovery::bit_count); ++integer) {
            BOOST_TEST_CONTEXT("modulus " << Modulus << ", integer " << integer) {
                const auto witness = witness_from_integer(cs, typename field::integral_type(integer), 1);
                BOOST_CHECK_EQUAL(rows_satisfied(cs, witness[0], witness, 0, cs.constraints.size()), integer < Modulus);
                BOOST_CHECK_EQUAL(cs.is_satisfied(witness[0], witness), integer < Modulus && (integer & 1) != 0);
                if (integer < Modulus && (integer & 1) != 0) {
                    std::vector<value> generated = {value(integer)};
                    recovery::append_witness(generated);
                    BOOST_CHECK(generated == witness);
                }
            }
        }
    }
}    // namespace

BOOST_AUTO_TEST_SUITE(r1cs_to_modified_sap_test_suite)

BOOST_AUTO_TEST_CASE(recovery_dimensions_and_source_prefix) {
    BOOST_CHECK_EQUAL(recovery_type::bit_count, 254);
    BOOST_CHECK_EQUAL(recovery_type::marker_count, 99);
    BOOST_CHECK_EQUAL(recovery_type::product_count, 153);
    BOOST_CHECK_EQUAL(recovery_type::witness_size, 506);
    BOOST_CHECK_EQUAL(recovery_type::constraint_count, 561);
    for (const std::size_t first_bit : {1, 3}) {
        const auto cs = recovery_system<field_type>(first_bit);
        BOOST_CHECK_EQUAL(cs.num_variables(), first_bit + 506);
        BOOST_CHECK_EQUAL(cs.num_constraints(), 562);
        BOOST_REQUIRE(cs.is_valid());
        const auto normalized = cs.normalized();
        for (const auto &u : {value_type(1), value_type(3), value_type(field_type::modulus - 2)}) {
            std::vector<value_type> witness(first_bit, value_type(13));
            witness[0] = u;
            const auto prefix = witness;
            recovery_type::append_witness(witness);
            auto expected = witness_from_integer(cs, u.to_integral(), first_bit);
            std::copy(prefix.begin(), prefix.end(), expected.begin());
            BOOST_CHECK(witness == expected);
            BOOST_CHECK_EQUAL_COLLECTIONS(witness.begin(), witness.begin() + first_bit, prefix.begin(), prefix.end());
            BOOST_CHECK_EQUAL(witness.size(), cs.num_variables());
            BOOST_CHECK(witness[first_bit].is_one());
            BOOST_CHECK(cs.is_satisfied(u, witness));
            BOOST_CHECK(normalized.is_satisfied(u, witness));
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_noncanonical_representatives_even_when_reconstruction_passes) {
    const auto cs = recovery_system<field_type>(1);
    for (const unsigned offset : {0, 1, 3}) {
        BOOST_TEST_CONTEXT("modulus plus " << offset) {
            const field_type::integral_type integer = field_type::modulus + offset;
            const auto witness = witness_from_integer(cs, integer, 1);
            BOOST_CHECK(witness[0] == value_type(offset));
            // Binding, Booleanity and field reconstruction all hold; comparison must reject the lift.
            BOOST_CHECK(rows_satisfied(cs, witness[0], witness, 0, recovery_type::bit_count + 2));
            BOOST_CHECK(!rows_satisfied(cs, witness[0], witness, recovery_type::bit_count + 2, cs.constraints.size()));
            BOOST_CHECK(!cs.is_satisfied(witness[0], witness));
            if (offset != 0) {
                BOOST_CHECK(witness[1].is_zero());
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(rejects_nonboolean_bits_even_when_reconstruction_and_comparison_pass) {
    const auto cs = recovery_system<field_type>(1);
    std::vector<value_type> witness(cs.num_variables(), value_type::zero());
    witness[0] = value_type(3);
    witness[1] = value_type(3);    // Represents 3 through a non-Boolean low bit.
    complete_comparison_witness(cs, witness, 1);
    BOOST_CHECK(rows_satisfied(cs, witness[0], witness, recovery_type::bit_count + 1, cs.constraints.size()));
    BOOST_CHECK(!rows_satisfied(cs, witness[0], witness, 1, recovery_type::bit_count + 1));
    BOOST_CHECK(!cs.is_satisfied(witness[0], witness));
}

BOOST_AUTO_TEST_CASE(rejects_changed_bits_markers_and_product_auxiliaries) {
    const auto cs = recovery_system<field_type>(1);
    const value_type u(field_type::modulus - 2);
    std::vector<value_type> witness = {u};
    recovery_type::append_witness(witness);
    BOOST_REQUIRE(cs.is_satisfied(u, witness));

    for (const std::size_t bit : {0, 1, 127, 253}) {
        auto changed = witness;
        changed[1 + bit] = value_type::one() - changed[1 + bit];
        complete_comparison_witness(cs, changed, 1);
        BOOST_CHECK(!cs.is_satisfied(u, changed));
    }
    for (std::size_t marker = 0; marker < recovery_type::marker_count; ++marker) {
        BOOST_TEST_CONTEXT("marker " << marker) {
            auto changed = witness;
            changed[1 + recovery_type::bit_count + marker] += value_type::one();
            complete_product_auxiliaries(cs, changed, 1);
            BOOST_CHECK(!cs.is_satisfied(u, changed));
        }
    }
    for (std::size_t auxiliary = 0; auxiliary < recovery_type::product_count; ++auxiliary) {
        BOOST_TEST_CONTEXT("product auxiliary " << auxiliary) {
            auto changed = witness;
            changed[1 + recovery_type::bit_count + recovery_type::marker_count + auxiliary] += value_type::one();
            BOOST_CHECK(!cs.is_satisfied(u, changed));
        }
    }
    auto changed = witness;
    changed[0] += value_type::one();
    BOOST_CHECK(!rows_satisfied(cs, u, changed, 0, 1));
    BOOST_CHECK(!cs.is_satisfied(u, changed));
    BOOST_CHECK(!cs.is_satisfied(value_type(1), witness));
}

BOOST_AUTO_TEST_CASE(invalid_inputs_fail_before_mutation) {
    for (const auto &u : {value_type(0), value_type(2), value_type(field_type::modulus - 1)}) {
        std::vector<value_type> witness = {u, value_type(13)};
        const auto original = witness;
        BOOST_CHECK_THROW(recovery_type::append_witness(witness), std::invalid_argument);
        BOOST_CHECK(witness == original);
    }
    std::vector<value_type> empty;
    BOOST_CHECK_THROW(recovery_type::append_witness(empty), std::invalid_argument);
    BOOST_CHECK(empty.empty());

    auto cs = recovery_system<field_type>(1);
    const auto original = cs.constraints;
    const auto max_size = std::numeric_limits<std::size_t>::max();
    for (const auto first_bit : {std::size_t(0), max_size, max_size - recovery_type::witness_size + 1,
                                 std::vector<value_type>().max_size() - recovery_type::witness_size + 1}) {
        BOOST_CHECK_THROW(recovery_type::append_constraints(cs.constraints, first_bit), std::invalid_argument);
        BOOST_CHECK(cs.constraints == original);
    }
}

BOOST_AUTO_TEST_CASE(exhaustive_small_field_bit_patterns) {
    test_all_bit_patterns<3>();
    test_all_bit_patterns<5>();
    test_all_bit_patterns<7>();
    test_all_bit_patterns<11>();
    test_all_bit_patterns<13>();
    test_all_bit_patterns<17>();
    test_all_bit_patterns<31>();
}

BOOST_AUTO_TEST_CASE(exhaustive_bits_and_marker_over_field_seven) {
    using field = test_prime_field<7>;
    using value = field::value_type;
    const auto cs = recovery_system<field>(1);
    // Exhaust field values for all three bits and the single marker; product auxiliaries are forced by the minus rows.
    for (unsigned b0 = 0; b0 < 7; ++b0) {
        for (unsigned b1 = 0; b1 < 7; ++b1) {
            for (unsigned b2 = 0; b2 < 7; ++b2) {
                for (unsigned marker = 0; marker < 7; ++marker) {
                    BOOST_TEST_CONTEXT("bits " << b0 << ',' << b1 << ',' << b2 << ", marker " << marker) {
                        std::vector<value> witness(cs.num_variables(), value::zero());
                        witness[0] = value(b0 + 2 * b1 + 4 * b2);
                        witness[1] = value(b0);
                        witness[2] = value(b1);
                        witness[3] = value(b2);
                        witness[4] = value(marker);
                        complete_product_auxiliaries(cs, witness, 1);
                        const bool valid =
                            b0 <= 1 && b1 <= 1 && b2 <= 1 && b0 + 2 * b1 + 4 * b2 < 7 && marker == b2 * b1;
                        BOOST_CHECK_EQUAL(rows_satisfied(cs, witness[0], witness, 0, cs.constraints.size()), valid);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
