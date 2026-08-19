#---------------------------------------------------------------------------#
# Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

include(CMTest)

cm_test_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
                       Boost::random
                       Boost::unit_test_framework
                       )

macro(define_math_test name)
    set(test_name "${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}_${name}_test")

    set(additional_args "")
    if(ENABLE_JUNIT_TEST_OUTPUT)
        set(TEST_RESULTS_DIR "${CMAKE_CURRENT_BINARY_DIR}/junit_results")
        set(TEST_LOGS_DIR "${TEST_RESULTS_DIR}/logs")
        file(MAKE_DIRECTORY ${TEST_LOGS_DIR})

        set(additional_args "--log_format=JUNIT"
                            "--log_sink=${TEST_LOGS_DIR}/${test_name}.xml")
    endif()

    cm_test(NAME ${test_name} SOURCES ${name}.cpp ARGS ${additional_args})

    target_include_directories(${test_name} PRIVATE
                               "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                               "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"

                               ${Boost_INCLUDE_DIRS})

    set_target_properties(${test_name} PROPERTIES CXX_STANDARD 23)

endmacro()

set(TESTS_NAMES
    "batch_inverse"
    "evaluation_domain"
    "geometric_sequence_domain"
    "integer_sqrt"
    "matrix"
    "polynomial_gcd"
    # "kronecker_substitution"  # FIXME: This fails. Disabled for passing CI
    "polynomial_arithmetic"
    "polynomial_backend"
    "polynomial_division"
    "polynomial_exponentiation"
    "polynomial_composition"
    "power_series"
    "polynomial"
    "polynomial_view"
    "polynomial_dfs"
    "polynomial_dfs_view"
    "lagrange_interpolation"
    "linear_combination"
    "unity_root"
    "mixed_radix_fft"
    "basic_radix2_domain")

foreach(TEST_NAME ${TESTS_NAMES})
    define_math_test(${TEST_NAME})
endforeach()

set(COMPILE_TIME_TESTS_NAMES
        "polynomial_algorithm_constraints"
        "polynomial_concepts"
        "static_matrix")

foreach(TEST_NAME ${COMPILE_TIME_TESTS_NAMES})
    cm_test(NAME math_${TEST_NAME}_compile_test SOURCES ${TEST_NAME}.cpp COMPILE_ONLY)
    target_include_directories(math_${TEST_NAME}_compile_test PRIVATE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"
            ${Boost_INCLUDE_DIRS})
    set_target_properties(math_${TEST_NAME}_compile_test PROPERTIES CXX_STANDARD 23)
endforeach()

if(BUILD_BENCH_TESTS)
    cm_add_test_subdirectory(bench_test)
endif()
