#---------------------------------------------------------------------------#
# Copyright (c) 2018-2021 Mikhail Komarov <nemo@nil.foundation>
# Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

include(CMTest)

add_custom_target(crypto3_benchmarks)

macro(define_benchmark benchmark)

    get_filename_component(name ${benchmark} NAME)
    string(REPLACE "/" "_" full_name ${benchmark}_benchmark)

    add_dependencies(crypto3_benchmarks ${full_name})

    cm_test(NAME ${full_name} SOURCES ${benchmark}.cpp)

    target_include_directories(
        ${full_name} PRIVATE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"
            ${Boost_INCLUDE_DIRS})

    target_link_libraries(${full_name}
        ${CMAKE_WORKSPACE_NAME}::random
        ${CMAKE_WORKSPACE_NAME}::math
        ${CMAKE_WORKSPACE_NAME}::algebra
        ${CMAKE_WORKSPACE_NAME}::multiprecision
        ${CMAKE_WORKSPACE_NAME}::benchmark_tools

        Boost::unit_test_framework
        Boost::timer
        Boost::random)

    if (TARGET ${CMAKE_WORKSPACE_NAME}::zk)
        target_link_libraries(${full_name} ${CMAKE_WORKSPACE_NAME}::zk)
    endif()

    set_target_properties(${full_name}
        PROPERTIES
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED TRUE)

    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(${full_name} PRIVATE "-fconstexpr-steps=2147483647")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${full_name} PRIVATE "-fconstexpr-ops-limit=4294967295")
    endif()

endmacro()

macro(define_standalone_benchmark benchmark)

    string(REPLACE "/" "_" full_name ${benchmark}_benchmark)

    add_dependencies(crypto3_benchmarks ${full_name})

    add_executable(${full_name} ${benchmark}.cpp)

    target_include_directories(
        ${full_name} PRIVATE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"
            ${Boost_INCLUDE_DIRS})

    target_link_libraries(${full_name}
        ${CMAKE_WORKSPACE_NAME}::random
        ${CMAKE_WORKSPACE_NAME}::math
        ${CMAKE_WORKSPACE_NAME}::algebra
        ${CMAKE_WORKSPACE_NAME}::multiprecision

        Boost::random)

    set_target_properties(${full_name}
        PROPERTIES
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED TRUE)

    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(${full_name} PRIVATE "-fconstexpr-steps=2147483647")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${full_name} PRIVATE "-fconstexpr-ops-limit=4294967295")
    endif()

endmacro()

set(BENCHMARK_NAMES
    "algebra/curves"
    "algebra/fields"
    "algebra/multiexp"

    "math/polynomial_dfs"

    "multiprecision/big_mod"

    "zk/lpc"
)

foreach(BENCHMARK_NAME ${BENCHMARK_NAMES})
    define_benchmark(${BENCHMARK_NAME})
endforeach()

set(STANDALONE_BENCHMARK_NAMES
    "algebra/bn254_fp12"
    "algebra/bn254_fp12_asm"
)

foreach(BENCHMARK_NAME ${STANDALONE_BENCHMARK_NAMES})
    define_standalone_benchmark(${BENCHMARK_NAME})
endforeach()

set(BN254_FP12_ASM_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/algebra/bn254_fp12_asm.cpp")
set(BN254_FP12_ASM_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/assembly/bn254_fp12_asm.s")

set(BN254_FP12_ASM_COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")
separate_arguments(BN254_FP12_ASM_FLAGS NATIVE_COMMAND "${BN254_FP12_ASM_COMPILE_FLAGS}")

string(TOUPPER "${CMAKE_BUILD_TYPE}" BN254_FP12_ASM_BUILD_TYPE)
if (BN254_FP12_ASM_BUILD_TYPE)
    set(BN254_FP12_ASM_CONFIG_FLAGS_RAW "${CMAKE_CXX_FLAGS_${BN254_FP12_ASM_BUILD_TYPE}}")
    separate_arguments(BN254_FP12_ASM_CONFIG_FLAGS NATIVE_COMMAND "${BN254_FP12_ASM_CONFIG_FLAGS_RAW}")
    list(APPEND BN254_FP12_ASM_FLAGS ${BN254_FP12_ASM_CONFIG_FLAGS})
endif()

list(APPEND BN254_FP12_ASM_FLAGS
    "-O3"
    "-g"
    "-DNDEBUG"
    "-std=gnu++23"
    "-S"
    "-fverbose-asm")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    list(APPEND BN254_FP12_ASM_FLAGS "-fconstexpr-steps=2147483647")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    list(APPEND BN254_FP12_ASM_FLAGS "-fconstexpr-ops-limit=4294967295")
endif()

if (CMAKE_TARGET_ARCHITECTURE MATCHES "x86|X86|amd64|AMD64")
    list(APPEND BN254_FP12_ASM_FLAGS "-masm=intel")
endif()

set(BN254_FP12_ASM_INCLUDE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_BINARY_DIR}/include"
    "${CMAKE_SOURCE_DIR}/libs/random/include"
    "${CMAKE_SOURCE_DIR}/libs/math/include"
    "${CMAKE_SOURCE_DIR}/libs/algebra/include"
    "${CMAKE_SOURCE_DIR}/libs/multiprecision/include"
    "${CMAKE_SOURCE_DIR}/libs/hash/include"
    "${CMAKE_SOURCE_DIR}/libs/marshalling/algebra/include"
    "${CMAKE_SOURCE_DIR}/libs/marshalling/multiprecision/include"
    "${CMAKE_SOURCE_DIR}/libs/marshalling/core/include"
    "${CMAKE_SOURCE_DIR}/libs/block/include"
    "${CMAKE_SOURCE_DIR}/libs/mac/include"
    "${CMAKE_SOURCE_DIR}/libs/parallelization-utils/include"
    "${CMAKE_SOURCE_DIR}/libs/benchmark_tools/include"
    ${Boost_INCLUDE_DIRS})

foreach(BN254_FP12_ASM_INCLUDE_DIR IN LISTS BN254_FP12_ASM_INCLUDE_DIRS)
    if (BN254_FP12_ASM_INCLUDE_DIR)
        list(APPEND BN254_FP12_ASM_INCLUDE_FLAGS "-I${BN254_FP12_ASM_INCLUDE_DIR}")
    endif()
endforeach()

add_custom_target(algebra_bn254_fp12_asm_assembly
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/assembly"
    COMMAND ${CMAKE_CXX_COMPILER}
            ${BN254_FP12_ASM_FLAGS}
            ${BN254_FP12_ASM_INCLUDE_FLAGS}
            -o "${BN254_FP12_ASM_OUTPUT}"
            "${BN254_FP12_ASM_SOURCE}"
    BYPRODUCTS "${BN254_FP12_ASM_OUTPUT}"
    DEPENDS "${BN254_FP12_ASM_SOURCE}"
    VERBATIM
    COMMAND_EXPAND_LISTS
    COMMENT "Generating BN254 fp12 assembly at ${BN254_FP12_ASM_OUTPUT}")
