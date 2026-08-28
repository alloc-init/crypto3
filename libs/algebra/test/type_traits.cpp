//---------------------------------------------------------------------------//
// Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
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

#define BOOST_TEST_MODULE type_traits_test

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/container/vector.hpp>

#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/algebra/curves/mnt4.hpp>
#include <nil/crypto3/algebra/curves/mnt6.hpp>
#include <nil/crypto3/algebra/curves/vesta.hpp>
#include <nil/crypto3/algebra/curves/pallas.hpp>
#include <nil/crypto3/algebra/curves/alt_bn128.hpp>
#include <nil/crypto3/algebra/curves/jubjub.hpp>
#include <nil/crypto3/algebra/curves/babyjubjub.hpp>
#include <nil/crypto3/algebra/curves/secp_k1.hpp>
#include <nil/crypto3/algebra/curves/secp_r1.hpp>
#include <nil/crypto3/algebra/curves/ed25519.hpp>

#include <nil/crypto3/algebra/fields/fp2.hpp>
#include <nil/crypto3/algebra/fields/fp3.hpp>
#include <nil/crypto3/algebra/fields/fp4.hpp>
#include <nil/crypto3/algebra/fields/fp6_3over2.hpp>
#include <nil/crypto3/algebra/fields/fp6_2over3.hpp>
#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>

#include <nil/crypto3/algebra/fields/babybear/base_field.hpp>
#include <nil/crypto3/algebra/fields/goldilocks.hpp>
#include <nil/crypto3/algebra/fields/koalabear.hpp>
#include <nil/crypto3/algebra/fields/mersenne31.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/mersenne31.hpp>

#include <nil/crypto3/algebra/type_traits.hpp>

using namespace nil::crypto3::algebra;

BOOST_AUTO_TEST_SUITE(type_traits_manual_tests)
/**/

template<typename value_type>
void test_field_value_types() {
    static_assert(FieldElementWithCoordinates<value_type>);

    static_assert(FieldValue<value_type>);
}

template<typename field_type>
void test_field_types() {
    static_assert(Field<field_type>);

    test_field_value_types<typename field_type::value_type>();
}

template<typename field_type>
void test_extended_field_types() {
    test_field_types<field_type>();

    static_assert(ExtendedField<field_type>);
    static_assert(ExtendedFieldValue<typename field_type::value_type>);

    test_field_value_types<typename field_type::value_type>();
}

template<typename curve_group_type>
void test_curve_group_types() {
    static_assert(CurveGroup<curve_group_type>);
    using value_type = typename curve_group_type::value_type;

    static_assert(CurveElement<value_type>);
}

template<typename curve_type>
void test_ordinary_curve_types() {
    test_field_types<typename curve_type::base_field_type>();

    test_field_types<typename curve_type::scalar_field_type>();

    test_curve_group_types<typename curve_type::template g1_type<>>();

    static_assert(Curve<curve_type>);
}

template<typename curve_type>
void test_pairing_friendly_curve_types() {
    test_ordinary_curve_types<curve_type>();

    static_assert(CurveWithG2<curve_type>);
    test_curve_group_types<typename curve_type::template g2_type<>>();

    using g2_base_field = typename curve_type::template g2_type<>::params_type::field_type;
    test_extended_field_types<g2_base_field>();

    static_assert(CurveWithTargetGroup<curve_type>);
    test_extended_field_types<typename curve_type::gt_type>();
}

BOOST_AUTO_TEST_CASE(pasta_type_traits) {
    test_ordinary_curve_types<curves::pallas>();
    test_ordinary_curve_types<curves::vesta>();
}

BOOST_AUTO_TEST_CASE(bls12_type_traits) {
    test_pairing_friendly_curve_types<curves::bls12<381>>();
    test_pairing_friendly_curve_types<curves::bls12<377>>();
}

BOOST_AUTO_TEST_CASE(mnt_type_traits) {
    test_pairing_friendly_curve_types<curves::mnt4<298>>();
    test_pairing_friendly_curve_types<curves::mnt6<298>>();
}

BOOST_AUTO_TEST_CASE(alt_bn128_type_traits) {
    test_pairing_friendly_curve_types<curves::alt_bn128<254>>();
    using base_field_type = curves::alt_bn128<254>::base_field_type;
    using fp2_field_type = fields::fp2<base_field_type>;
    using fp6_field_type = fields::fp6_3over2<base_field_type>;
    using fp12_field_type = fields::fp12_2over3over2<base_field_type>;
    static_assert(FieldValue<typename base_field_type::value_type>);
    static_assert(FieldValue<typename fp2_field_type::value_type>);
    static_assert(FieldValue<typename fp6_field_type::value_type>);
    static_assert(FieldValue<typename fp12_field_type::value_type>);
    static_assert(!FieldValue<base_field_type>);
}

BOOST_AUTO_TEST_CASE(jubjub_type_traits) {
    test_ordinary_curve_types<curves::jubjub>();
}

BOOST_AUTO_TEST_CASE(babyjubjub_type_traits) {
    test_ordinary_curve_types<curves::babyjubjub>();
}

BOOST_AUTO_TEST_CASE(goldilocks_field_type_traits) {
    test_field_types<fields::goldilocks>();
    test_field_types<fields::goldilocks_fp2>();
}

BOOST_AUTO_TEST_CASE(mersenne31_field_type_traits) {
    test_field_types<fields::mersenne31>();

    using field_type = fields::mersenne31;
    using value_type = field_type::value_type;
    using params_type = fields::arithmetic_params<field_type>;

    static_assert(field_type::modulus == 0x7fffffffu);
    static_assert(params_type::s == 1);
    static_assert(value_type(params_type::root_of_unity) == -value_type::one());
}

BOOST_AUTO_TEST_CASE(koalabear_field_type_traits) {
    test_field_types<fields::koalabear>();
}

BOOST_AUTO_TEST_CASE(babybear_field_type_traits) {
    test_field_types<fields::babybear>();
    test_field_types<fields::babybear_fp4>();
}

BOOST_AUTO_TEST_CASE(secp_type_traits) {
    test_ordinary_curve_types<curves::secp160r1>();
    test_ordinary_curve_types<curves::secp192r1>();
    test_ordinary_curve_types<curves::secp224r1>();
    test_ordinary_curve_types<curves::secp256r1>();
    test_ordinary_curve_types<curves::secp384r1>();
    test_ordinary_curve_types<curves::secp521r1>();

    test_ordinary_curve_types<curves::secp160k1>();
    test_ordinary_curve_types<curves::secp192k1>();
    test_ordinary_curve_types<curves::secp224k1>();
    test_ordinary_curve_types<curves::secp256k1>();
}

BOOST_AUTO_TEST_CASE(ed25519_type_traits) {
    test_ordinary_curve_types<curves::ed25519>();
}

template<typename T>
concept HasSqrt = requires(const T &value) {
    { value.sqrt() } -> std::same_as<T>;
};

#define FIELD_HAS_SQRT(field) (HasSqrt<typename field::value_type>)

BOOST_AUTO_TEST_CASE(test_extended_fields_sqrt_trait) {

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::alt_bn128_254::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::alt_bn128_254::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::alt_bn128_254::template g1_type<>::field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::alt_bn128_254::template g2_type<>::field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::alt_bn128_254::gt_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_381::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_381::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_381::template g1_type<>::field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_381::template g2_type<>::field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_381::gt_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_377::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_377::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_377::template g1_type<>::field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_377::template g2_type<>::field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::bls12_377::gt_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::mnt4_298::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::mnt4_298::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::mnt4_298::template g1_type<>::field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::mnt4_298::template g2_type<>::field_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::mnt6_298::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::mnt6_298::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::mnt6_298::template g1_type<>::field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::mnt6_298::template g2_type<>::field_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::pallas::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::pallas::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::pallas::template g1_type<>::field_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::vesta::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::vesta::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::vesta::template g1_type<>::field_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::jubjub::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::jubjub::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::jubjub::template g1_type<>::field_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::babyjubjub::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::babyjubjub::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::babyjubjub::template g1_type<>::field_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(curves::ed25519::base_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::ed25519::scalar_field_type));
    BOOST_ASSERT(FIELD_HAS_SQRT(curves::ed25519::template g1_type<>::field_type));

    BOOST_ASSERT(FIELD_HAS_SQRT(fields::goldilocks));
    BOOST_ASSERT(FIELD_HAS_SQRT(fields::mersenne31));
    BOOST_ASSERT(FIELD_HAS_SQRT(fields::koalabear));
    BOOST_ASSERT(FIELD_HAS_SQRT(fields::babybear));
}

BOOST_AUTO_TEST_CASE(test_extended_fields_trait) {

    BOOST_ASSERT(!ExtendedFieldValue<curves::alt_bn128_254::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::alt_bn128_254::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::alt_bn128_254::template g1_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::alt_bn128_254::template g2_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::alt_bn128_254::gt_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::bls12_381::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::bls12_381::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::bls12_381::template g1_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::bls12_381::template g2_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::bls12_381::gt_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::bls12_377::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::bls12_377::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::bls12_377::template g1_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::bls12_377::template g2_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::bls12_377::gt_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::mnt4_298::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::mnt4_298::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::mnt4_298::template g1_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::mnt4_298::template g2_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::mnt4_298::gt_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::mnt6_298::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::mnt6_298::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::mnt6_298::template g1_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::mnt6_298::template g2_type<>::field_type::value_type>);
    BOOST_ASSERT(ExtendedFieldValue<curves::mnt6_298::gt_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::pallas::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::pallas::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::pallas::template g1_type<>::field_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::vesta::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::vesta::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::vesta::template g1_type<>::field_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::jubjub::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::jubjub::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::jubjub::template g1_type<>::field_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::babyjubjub::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::babyjubjub::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::babyjubjub::template g1_type<>::field_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<curves::ed25519::base_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::ed25519::scalar_field_type::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<curves::ed25519::template g1_type<>::field_type::value_type>);

    BOOST_ASSERT(!ExtendedFieldValue<fields::goldilocks::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<fields::mersenne31::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<fields::koalabear::value_type>);
    BOOST_ASSERT(!ExtendedFieldValue<fields::babybear::value_type>);
}

BOOST_AUTO_TEST_SUITE_END()
