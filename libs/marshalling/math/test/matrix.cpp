#define BOOST_TEST_MODULE crypto3_marshalling_matrix_test

#include <cstdint>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/marshalling/math/types/matrix.hpp>
#include <nil/crypto3/math/matrix/compressed.hpp>
#include <nil/marshalling/endianness.hpp>
#include <nil/marshalling/status_type.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;
    namespace math = nil::crypto3::math;
    namespace types = nil::crypto3::marshalling::types;

    using fp12_type = fields::fp12_2over3over2<fields::alt_bn128<254>>::value_type;
    using matrix_type = math::compressed_matrix<fp12_type>;

    template<typename Endianness>
    void check_round_trip() {
        matrix_type input(3, 5);
        input(0, 4) = fp12_type::one();
        input(2, 1) = fp12_type::one() + fp12_type::one();

        auto filled = types::fill_compressed_matrix<Endianness>(input);
        std::vector<std::uint8_t> bytes(filled.length());
        auto write_iter = bytes.begin();
        BOOST_REQUIRE(filled.write(write_iter, bytes.size()) == nil::marshalling::status_type::success);

        using filled_type = typename types::compressed_matrix<
            nil::marshalling::field_type<Endianness>, matrix_type>::type;
        filled_type decoded;
        auto read_iter = bytes.begin();
        BOOST_REQUIRE(decoded.read(read_iter, bytes.size()) == nil::marshalling::status_type::success);
        BOOST_REQUIRE(read_iter == bytes.end());

        const auto output = types::make_compressed_matrix<Endianness, matrix_type>(decoded);
        BOOST_CHECK_EQUAL(output.rows(), 3);
        BOOST_CHECK_EQUAL(output.columns(), 5);
        BOOST_REQUIRE(math::find_element(output, 0, 4) != nullptr);
        BOOST_REQUIRE(math::find_element(output, 2, 1) != nullptr);
        BOOST_CHECK(*math::find_element(output, 0, 4) == fp12_type::one());
        BOOST_CHECK(*math::find_element(output, 2, 1) == fp12_type::one() + fp12_type::one());
        BOOST_CHECK(math::find_element(output, 1, 1) == nullptr);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(compressed_fp12_matrix_big_endian) {
    check_round_trip<nil::marshalling::option::big_endian>();
}

BOOST_AUTO_TEST_CASE(compressed_fp12_matrix_little_endian) {
    check_round_trip<nil::marshalling::option::little_endian>();
}
