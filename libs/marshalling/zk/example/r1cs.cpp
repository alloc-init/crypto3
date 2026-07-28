#include <cstdint>
#include <iostream>
#include <vector>

#include <nil/marshalling/field_type.hpp>
#include <nil/marshalling/options.hpp>
#include <nil/marshalling/status_type.hpp>

#include <nil/crypto3/algebra/curves/bls12.hpp>
#include <nil/crypto3/marshalling/zk/types/r1cs_gg_ppzksnark/r1cs.hpp>

int main() {
    using field_type = nil::crypto3::algebra::curves::bls12<381>::scalar_field_type;
    using constraint_system_type = nil::crypto3::zk::snark::r1cs_constraint_system<field_type>;
    using endianness = nil::marshalling::option::big_endian;
    using type_base = nil::marshalling::field_type<endianness>;
    using marshalled_type =
        nil::crypto3::marshalling::types::r1cs_constraint_system<type_base, constraint_system_type>;

    using variable_type = nil::crypto3::zk::snark::linear_variable<field_type>;
    using linear_combination_type = nil::crypto3::zk::snark::linear_combination<variable_type>;

    constraint_system_type constraint_system;
    constraint_system.primary_input_size = 2;
    constraint_system.auxiliary_input_size = 1;

    linear_combination_type a, b, c;
    a.add_term(variable_type(1), field_type::value_type::one());
    b.add_term(variable_type(2), field_type::value_type::one());
    c.add_term(variable_type(3), field_type::value_type::one());
    constraint_system.add_constraint({a, b, c});

    auto filled = nil::crypto3::marshalling::types::fill_r1cs_constraint_system<
        constraint_system_type, endianness>(constraint_system);

    std::vector<std::uint8_t> bytes(filled.length());
    auto write_iter = bytes.begin();
    auto status = filled.write(write_iter, bytes.size());
    if (status != nil::marshalling::status_type::success) {
        std::cerr << "Failed to serialize R1CS constraint system" << std::endl;
        return 1;
    }

    marshalled_type decoded;
    auto read_iter = bytes.begin();
    status = decoded.read(read_iter, bytes.size());
    if (status != nil::marshalling::status_type::success) {
        std::cerr << "Failed to deserialize R1CS constraint system" << std::endl;
        return 1;
    }

    const auto restored = nil::crypto3::marshalling::types::make_r1cs_constraint_system<
        constraint_system_type, endianness>(decoded);
    if (restored != constraint_system) {
        std::cerr << "R1CS round trip changed the constraint system" << std::endl;
        return 1;
    }

    std::cout << "Serialized " << restored.num_constraints() << " R1CS constraints into "
              << bytes.size() << " bytes" << std::endl;
    return 0;
}
