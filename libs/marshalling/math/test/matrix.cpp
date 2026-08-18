#define BOOST_TEST_MODULE crypto3_marshalling_matrix_test

#include <cstdint>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/fp12_2over3over2.hpp>
#include <nil/crypto3/marshalling/math/types/matrix.hpp>
#include <nil/crypto3/math/matrix/compressed.hpp>
#include <nil/crypto3/math/matrix/regular.hpp>
#include <nil/marshalling/endianness.hpp>
#include <nil/marshalling/status_type.hpp>

namespace {
    namespace fields = nil::crypto3::algebra::fields;
    namespace math = nil::crypto3::math;
    namespace types = nil::crypto3::marshalling::types;

    using fp12_type = fields::fp12_2over3over2<fields::alt_bn128<254>>::value_type;
    using compressed_matrix_type = math::compressed_matrix<fp12_type>;
    using compressed_vector_type = math::compressed_vector<fp12_type>;
    using regular_matrix_type = math::regular_matrix<fp12_type>;
    using regular_vector_type = math::regular_vector<fp12_type>;

    template<typename FilledType>
    FilledType encode_and_decode(FilledType filled) {
        std::vector<std::uint8_t> bytes(filled.length());
        auto write_iter = bytes.begin();
        BOOST_REQUIRE(filled.write(write_iter, bytes.size()) == nil::marshalling::status_type::success);

        FilledType decoded;
        auto read_iter = bytes.begin();
        BOOST_REQUIRE(decoded.read(read_iter, bytes.size()) == nil::marshalling::status_type::success);
        BOOST_REQUIRE(read_iter == bytes.end());
        return decoded;
    }

    template<typename Endianness>
    void check_compressed_matrix_round_trip() {
        compressed_matrix_type input(3, 5);
        input(0, 4) = fp12_type::one();
        input(2, 1) = fp12_type::one() + fp12_type::one();

        auto filled = types::fill_compressed_matrix<Endianness>(input);
        const auto decoded = encode_and_decode(std::move(filled));
        const auto output = types::make_compressed_matrix<Endianness, compressed_matrix_type>(decoded);
        BOOST_CHECK_EQUAL(output.rows(), 3);
        BOOST_CHECK_EQUAL(output.columns(), 5);
        BOOST_REQUIRE(math::find_element(output, 0, 4) != nullptr);
        BOOST_REQUIRE(math::find_element(output, 2, 1) != nullptr);
        BOOST_CHECK(*math::find_element(output, 0, 4) == fp12_type::one());
        BOOST_CHECK(*math::find_element(output, 2, 1) == fp12_type::one() + fp12_type::one());
        BOOST_CHECK(math::find_element(output, 1, 1) == nullptr);
    }

    template<typename Endianness>
    void check_regular_matrix_round_trip() {
        regular_matrix_type input(2, 3);
        fp12_type value = fp12_type::one();
        for (std::size_t row = 0; row < input.rows(); ++row) {
            for (std::size_t column = 0; column < input.columns(); ++column) {
                input(row, column) = value;
                value += fp12_type::one();
            }
        }

        const auto decoded = encode_and_decode(types::fill_regular_matrix<Endianness>(input));
        const auto output = types::make_regular_matrix<Endianness, regular_matrix_type>(decoded);
        BOOST_CHECK_EQUAL(output.rows(), input.rows());
        BOOST_CHECK_EQUAL(output.columns(), input.columns());
        for (std::size_t row = 0; row < input.rows(); ++row) {
            for (std::size_t column = 0; column < input.columns(); ++column) {
                BOOST_CHECK(output(row, column) == input(row, column));
            }
        }
    }

    template<typename Endianness>
    void check_regular_vector_round_trip() {
        regular_vector_type input(3);
        input[0] = fp12_type::one();
        input[1] = fp12_type::one() + fp12_type::one();
        input[2] = fp12_type::one() + fp12_type::one() + fp12_type::one();

        const auto decoded = encode_and_decode(types::fill_regular_vector<Endianness>(input));
        const auto output = types::make_regular_vector<Endianness, regular_vector_type>(decoded);
        BOOST_REQUIRE_EQUAL(output.size(), input.size());
        for (std::size_t index = 0; index < input.size(); ++index) {
            BOOST_CHECK(output[index] == input[index]);
        }
    }

    template<typename Endianness>
    void check_compressed_vector_round_trip() {
        compressed_vector_type input(5);
        input[1] = fp12_type::one();
        input[4] = fp12_type::one() + fp12_type::one();

        const auto decoded = encode_and_decode(types::fill_compressed_vector<Endianness>(input));
        const auto output = types::make_compressed_vector<Endianness, compressed_vector_type>(decoded);
        BOOST_CHECK_EQUAL(output.size(), input.size());
        BOOST_REQUIRE(math::find_element(output, 1) != nullptr);
        BOOST_REQUIRE(math::find_element(output, 4) != nullptr);
        BOOST_CHECK(*math::find_element(output, 1) == input[1]);
        BOOST_CHECK(*math::find_element(output, 4) == input[4]);
        BOOST_CHECK(math::find_element(output, 2) == nullptr);
    }
}    // namespace

BOOST_AUTO_TEST_CASE(compressed_fp12_matrix_big_endian) {
    check_compressed_matrix_round_trip<nil::marshalling::option::big_endian>();
}

BOOST_AUTO_TEST_CASE(compressed_fp12_matrix_little_endian) {
    check_compressed_matrix_round_trip<nil::marshalling::option::little_endian>();
}

BOOST_AUTO_TEST_CASE(other_fp12_matrix_types_big_endian) {
    check_regular_matrix_round_trip<nil::marshalling::option::big_endian>();
    check_regular_vector_round_trip<nil::marshalling::option::big_endian>();
    check_compressed_vector_round_trip<nil::marshalling::option::big_endian>();
}

BOOST_AUTO_TEST_CASE(other_fp12_matrix_types_little_endian) {
    check_regular_matrix_round_trip<nil::marshalling::option::little_endian>();
    check_regular_vector_round_trip<nil::marshalling::option::little_endian>();
    check_compressed_vector_round_trip<nil::marshalling::option::little_endian>();
}
