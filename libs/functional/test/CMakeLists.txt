#---------------------------------------------------------------------------#
# Copyright (c) 2018-2024 Mikhail Komarov <nemo@allocin.it>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

include(CMTest)

cm_test_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
        ${${CURRENT_PROJECT_NAME}_INTERFACE_LIBRARIES}
        Boost::unit_test_framework)

macro(define_functional_test test)
    get_filename_component(name ${test} NAME)
    string(REPLACE "/" "_" test_name ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}_${test}_test)

    cm_test(NAME ${test_name} SOURCES ${test}.cpp)

    target_include_directories(${test_name} PRIVATE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"

            ${Boost_INCLUDE_DIRS})

    set_target_properties(${test_name} PROPERTIES
            CXX_STANDARD 20
            CXX_STANDARD_REQUIRED TRUE)

    target_compile_options(${test_name} PRIVATE "-ftemplate-backtrace-limit=0")

    string(CONCAT TEST_DATA ${CMAKE_CURRENT_SOURCE_DIR} "/data/" "${name}" ".json")
    target_compile_definitions(${test_name} PRIVATE TEST_DATA="${TEST_DATA}")

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(${test_name} PRIVATE "-fconstexpr-steps=2147483647")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${test_name} PRIVATE "-fconstexpr-ops-limit=4294967295")
    endif()

endmacro()

set(TESTS_NAMES
#        "abe/dippe"
#        "abe/fame"
#        "abe/gpsw"
#        "abe/policy"
#
#        "ipe/damgard"
#        "ipe/damgard_dec_multi"
#        "ipe/damgard_multi"
#        "ipe/dmcfe"
#        "ipe/fh_multi_ipe"
#        "ipe/fhipe"
#        "ipe/lwe_fs"
#        "ipe/paillier"
#
#        "ipe/ddh"
#        "ipe/ddh_multi"
#        "ipe/lwe"
#        "ipe/ring_lwe"
#
#        "internal/big"
#        "internal/dlog"
#        "internal/keygen"
#        "internal/prime"
#        "internal/str"
#
#        "quadratic/sgp"
#
#        "sample/normal"
#        "sample/normal_cdt"
#        "sample/normal_cumulative"
#        "sample/normal_double"
#        "sample/normal_double_constant"
#        "sample/normal_negative"
)

foreach(TEST_NAME ${TESTS_NAMES})
    define_functional_test(${TEST_NAME})
endforeach()