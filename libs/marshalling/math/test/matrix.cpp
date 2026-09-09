#define BOOST_TEST_MODULE crypto3_marshalling_matrix_test

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
#include <list>
#include <tuple>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/fields/alt_bn128/base_field.hpp>
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
    using base_field_type = fields::alt_bn128_base_field<254>::value_type;
    using compressed_matrix_type = math::compressed_matrix<fp12_type>;
    using compressed_vector_type = math::compressed_vector<fp12_type>;
    using base_regular_matrix_type = math::regular_matrix<base_field_type>;
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

    template<typename Endianness>
    void check_incremental_regular_matrix_encoding(const regular_matrix_type &matrix) {
        const auto filled = types::fill_regular_matrix<Endianness>(matrix);
        std::vector<std::uint8_t> expected(filled.length());
        auto expected_output = expected.begin();
        BOOST_REQUIRE(filled.write(expected_output, expected.size()) == nil::marshalling::status_type::success);

        std::vector<std::uint8_t> incremental;
        incremental.reserve(expected.size());
        auto incremental_output = std::back_inserter(incremental);
        BOOST_REQUIRE((types::write_regular_matrix<Endianness>(matrix, incremental_output, expected.size()) ==
                       nil::marshalling::status_type::success));
        BOOST_CHECK(incremental == expected);
    }

    template<typename Endianness, typename MatrixType>
    std::vector<std::uint8_t> encode_regular_matrix(const MatrixType &matrix) {
        const auto filled = types::fill_regular_matrix<Endianness>(matrix);
        std::vector<std::uint8_t> bytes(filled.length());
        auto output = bytes.begin();
        BOOST_REQUIRE(filled.write(output, bytes.size()) == nil::marshalling::status_type::success);
        return bytes;
    }

    class counting_byte_output_iterator {
    public:
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::output_iterator_tag;
        using pointer = void;
        using reference = void;
        using value_type = void;

        explicit counting_byte_output_iterator(std::size_t &count) : count_(&count) {
        }

        counting_byte_output_iterator &operator*() {
            return *this;
        }

        counting_byte_output_iterator &operator=(std::uint8_t) {
            ++*count_;
            return *this;
        }

        counting_byte_output_iterator &operator++() {
            return *this;
        }

        counting_byte_output_iterator operator++(int) {
            return *this;
        }

    private:
        std::size_t *count_;
    };

    class generated_base_field_matrix {
    public:
        using value_type = base_field_type;

        generated_base_field_matrix(std::size_t rows, std::size_t columns) : rows_(rows), columns_(columns) {
        }

        std::size_t rows() const {
            return rows_;
        }

        std::size_t columns() const {
            return columns_;
        }

        const value_type &operator()(std::size_t, std::size_t) const {
            ++element_accesses_;
            return value_;
        }

        std::size_t element_accesses() const {
            return element_accesses_;
        }

    private:
        std::size_t rows_;
        std::size_t columns_;
        value_type value_ = value_type::one();
        mutable std::size_t element_accesses_ = 0;
    };
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

BOOST_AUTO_TEST_CASE(regular_matrix_encoded_size_is_checked_without_materializing_elements) {
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using marshalling_type = types::regular_matrix<type_base, regular_matrix_type>;
    constexpr std::size_t header_length =
        2 * marshalling_type::index_type::max_length() + marshalling_type::element_count_type::max_length();

    std::size_t encoded_size = 0;
    BOOST_REQUIRE((types::regular_matrix_encoded_size<endianness, regular_matrix_type>(0, 7, encoded_size) ==
                   nil::marshalling::status_type::success));
    BOOST_CHECK_EQUAL(encoded_size, header_length);

    regular_matrix_type rectangular_matrix(2, 3);
    const auto filled_matrix = types::fill_regular_matrix<endianness>(rectangular_matrix);
    BOOST_REQUIRE((types::regular_matrix_encoded_size<endianness, regular_matrix_type>(2, 3, encoded_size) ==
                   nil::marshalling::status_type::success));
    BOOST_CHECK_EQUAL(encoded_size, filled_matrix.length());

    constexpr std::size_t maximum_size = std::numeric_limits<std::size_t>::max();
    encoded_size = 17;
    BOOST_CHECK((types::regular_matrix_encoded_size<endianness, regular_matrix_type>(maximum_size, 2, encoded_size) ==
                 nil::marshalling::status_type::invalid_msg_data));
    BOOST_CHECK_EQUAL(encoded_size, 17);

    constexpr std::size_t element_length = marshalling_type::value_type::max_length();
    const std::size_t element_count_with_oversized_encoding = maximum_size / element_length + 1;
    BOOST_CHECK((types::regular_matrix_encoded_size<endianness, regular_matrix_type>(
                     element_count_with_oversized_encoding, 1, encoded_size) ==
                 nil::marshalling::status_type::invalid_msg_data));
    BOOST_CHECK_EQUAL(encoded_size, 17);
}

BOOST_AUTO_TEST_CASE(incremental_regular_matrix_encoding_matches_the_existing_wire_format) {
    regular_matrix_type rectangular(2, 3);
    fp12_type value = fp12_type::one();
    for (std::size_t row = 0; row < rectangular.rows(); ++row) {
        for (std::size_t column = 0; column < rectangular.columns(); ++column) {
            rectangular(row, column) = value;
            value += fp12_type::one();
        }
    }
    check_incremental_regular_matrix_encoding<nil::marshalling::option::big_endian>(rectangular);
    check_incremental_regular_matrix_encoding<nil::marshalling::option::little_endian>(rectangular);

    const regular_matrix_type empty(0, 7);
    check_incremental_regular_matrix_encoding<nil::marshalling::option::big_endian>(empty);
    check_incremental_regular_matrix_encoding<nil::marshalling::option::little_endian>(empty);
}

BOOST_AUTO_TEST_CASE(incremental_regular_matrix_encoding_streams_without_matrix_or_output_materialization) {
    using endianness = nil::marshalling::option::big_endian;
    const generated_base_field_matrix generated(256, 256);
    std::size_t encoded_size = 0;
    BOOST_REQUIRE((types::regular_matrix_encoded_size<endianness, generated_base_field_matrix>(
                       generated.rows(), generated.columns(), encoded_size) == nil::marshalling::status_type::success));

    std::size_t written_bytes = 0;
    counting_byte_output_iterator output(written_bytes);
    BOOST_REQUIRE((types::write_regular_matrix<endianness>(generated, output, encoded_size) ==
                   nil::marshalling::status_type::success));
    BOOST_CHECK_EQUAL(generated.element_accesses(), generated.rows() * generated.columns());
    BOOST_CHECK_EQUAL(written_bytes, encoded_size);
}

BOOST_AUTO_TEST_CASE(incremental_regular_matrix_encoding_reports_size_errors_before_writing) {
    using endianness = nil::marshalling::option::big_endian;
    const generated_base_field_matrix one_element(1, 1);
    std::size_t encoded_size = 0;
    BOOST_REQUIRE((types::regular_matrix_encoded_size<endianness, generated_base_field_matrix>(1, 1, encoded_size) ==
                   nil::marshalling::status_type::success));

    std::size_t written_bytes = 0;
    counting_byte_output_iterator output(written_bytes);
    BOOST_CHECK((types::write_regular_matrix<endianness>(one_element, output, encoded_size - 1) ==
                 nil::marshalling::status_type::buffer_overflow));
    BOOST_CHECK_EQUAL(written_bytes, 0);
    BOOST_CHECK_EQUAL(one_element.element_accesses(), 0);

    const generated_base_field_matrix overflowing(std::numeric_limits<std::size_t>::max(), 2);
    BOOST_CHECK(
        (types::write_regular_matrix<endianness>(overflowing, output, std::numeric_limits<std::size_t>::max()) ==
         nil::marshalling::status_type::invalid_msg_data));
    BOOST_CHECK_EQUAL(written_bytes, 0);
    BOOST_CHECK_EQUAL(overflowing.element_accesses(), 0);
}

BOOST_AUTO_TEST_CASE(incremental_regular_matrix_decoding_visits_values_in_row_major_order) {
    using endianness = nil::marshalling::option::big_endian;
    regular_matrix_type matrix(2, 3);
    fp12_type value = fp12_type::one();
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            matrix(row, column) = value;
            value += fp12_type::one();
        }
    }

    const std::vector<std::uint8_t> bytes = encode_regular_matrix<endianness>(matrix);
    const std::list<std::uint8_t> sequential_bytes(bytes.begin(), bytes.end());
    auto input = sequential_bytes.begin();
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::vector<std::tuple<std::size_t, std::size_t, fp12_type>> visited;
    BOOST_REQUIRE((types::read_regular_matrix<endianness, regular_matrix_type>(
                       input, sequential_bytes.size(), rows, columns,
                       [&visited](std::size_t row, std::size_t column, const fp12_type &element) {
                           visited.emplace_back(row, column, element);
                       }) == nil::marshalling::status_type::success));

    BOOST_CHECK(input == sequential_bytes.end());
    BOOST_CHECK_EQUAL(rows, matrix.rows());
    BOOST_CHECK_EQUAL(columns, matrix.columns());
    BOOST_REQUIRE_EQUAL(visited.size(), matrix.rows() * matrix.columns());
    std::size_t index = 0;
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t column = 0; column < matrix.columns(); ++column, ++index) {
            BOOST_CHECK_EQUAL(std::get<0>(visited[index]), row);
            BOOST_CHECK_EQUAL(std::get<1>(visited[index]), column);
            BOOST_CHECK(std::get<2>(visited[index]) == matrix(row, column));
        }
    }
}

BOOST_AUTO_TEST_CASE(incremental_regular_matrix_decoding_accepts_an_empty_matrix) {
    using endianness = nil::marshalling::option::little_endian;
    const regular_matrix_type matrix(0, 7);
    const std::vector<std::uint8_t> bytes = encode_regular_matrix<endianness>(matrix);
    auto input = bytes.begin();
    std::size_t rows = 1;
    std::size_t columns = 1;
    std::size_t visits = 0;
    BOOST_REQUIRE((types::read_regular_matrix<endianness, regular_matrix_type>(
                       input, bytes.size(), rows, columns, [&visits](std::size_t, std::size_t, const fp12_type &) {
                           ++visits;
                       }) == nil::marshalling::status_type::success));
    BOOST_CHECK(input == bytes.end());
    BOOST_CHECK_EQUAL(rows, 0);
    BOOST_CHECK_EQUAL(columns, 7);
    BOOST_CHECK_EQUAL(visits, 0);
}

BOOST_AUTO_TEST_CASE(incremental_regular_matrix_decoding_rejects_truncation_and_an_incorrect_element_count) {
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using marshalling_type = types::regular_matrix<type_base, regular_matrix_type>;
    const regular_matrix_type matrix(2, 3);
    std::vector<std::uint8_t> bytes = encode_regular_matrix<endianness>(matrix);

    auto truncated_input = bytes.begin();
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::size_t visits = 0;
    BOOST_CHECK((types::read_regular_matrix<endianness, regular_matrix_type>(
                     truncated_input, bytes.size() - 1, rows, columns,
                     [&visits](std::size_t, std::size_t, const fp12_type &) { ++visits; }) ==
                 nil::marshalling::status_type::not_enough_data));
    BOOST_CHECK_EQUAL(visits, 0);

    const typename marshalling_type::element_count_type incorrect_element_count(5);
    auto count_output = bytes.begin() + 2 * marshalling_type::index_type::max_length();
    BOOST_REQUIRE(incorrect_element_count.write(count_output, incorrect_element_count.length()) ==
                  nil::marshalling::status_type::success);
    auto incorrect_count_input = bytes.begin();
    BOOST_CHECK((types::read_regular_matrix<endianness, regular_matrix_type>(
                     incorrect_count_input, bytes.size(), rows, columns,
                     [&visits](std::size_t, std::size_t, const fp12_type &) { ++visits; }) ==
                 nil::marshalling::status_type::invalid_msg_data));
    BOOST_CHECK_EQUAL(visits, 0);
}

BOOST_AUTO_TEST_CASE(incremental_regular_matrix_decoding_rejects_dimension_overflow) {
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using marshalling_type = types::regular_matrix<type_base, regular_matrix_type>;
    constexpr std::size_t header_length =
        2 * marshalling_type::index_type::max_length() + marshalling_type::element_count_type::max_length();
    std::vector<std::uint8_t> bytes(header_length);
    auto output = bytes.begin();
    const typename marshalling_type::index_type rows_field(std::numeric_limits<std::size_t>::max());
    const typename marshalling_type::index_type columns_field(2);
    const typename marshalling_type::element_count_type element_count_field(0);
    BOOST_REQUIRE(rows_field.write(output, rows_field.length()) == nil::marshalling::status_type::success);
    BOOST_REQUIRE(columns_field.write(output, columns_field.length()) == nil::marshalling::status_type::success);
    BOOST_REQUIRE(element_count_field.write(output, element_count_field.length()) ==
                  nil::marshalling::status_type::success);

    auto input = bytes.begin();
    std::size_t rows = 0;
    std::size_t columns = 0;
    BOOST_CHECK((types::read_regular_matrix<endianness, regular_matrix_type>(
                     input, bytes.size(), rows, columns, [](std::size_t, std::size_t, const fp12_type &) { }) ==
                 nil::marshalling::status_type::invalid_msg_data));
}

BOOST_AUTO_TEST_CASE(incremental_regular_matrix_decoding_rejects_noncanonical_field_elements) {
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using marshalling_type = types::regular_matrix<type_base, base_regular_matrix_type>;
    const base_regular_matrix_type matrix(1, 1);
    std::vector<std::uint8_t> bytes = encode_regular_matrix<endianness>(matrix);
    constexpr std::size_t header_length =
        2 * marshalling_type::index_type::max_length() + marshalling_type::element_count_type::max_length();
    std::fill(bytes.begin() + header_length, bytes.end(), 0xff);

    auto input = bytes.begin();
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::size_t visits = 0;
    BOOST_CHECK((types::read_regular_matrix<endianness, base_regular_matrix_type>(
                     input, bytes.size(), rows, columns, [&visits](std::size_t, std::size_t, const base_field_type &) {
                         ++visits;
                     }) == nil::marshalling::status_type::invalid_msg_data));
    BOOST_CHECK_EQUAL(visits, 0);
}
