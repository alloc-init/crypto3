#---------------------------------------------------------------------------#
# Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

include(CMTest)

cm_test_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
                       ${CMAKE_WORKSPACE_NAME}::multiprecision

                       Boost::unit_test_framework)

macro(define_codec_test name)
    set(TEST_TARGET "codec_${name}_test")

    cm_test(NAME ${TEST_TARGET} SOURCES ${name}.cpp)

    target_include_directories(${TEST_TARGET} PRIVATE
                               "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                               "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"

                               ${Boost_INCLUDE_DIRS})

    if(NOT CMAKE_CXX_STANDARD)
        set_target_properties(${TEST_TARGET} PROPERTIES
                              CXX_STANDARD 14
                              CXX_STANDARD_REQUIRED TRUE)
    endif()

    target_compile_definitions(${TEST_TARGET} PRIVATE TEST_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data")

endmacro()

set(TESTS_NAMES "hex" "base")

foreach(TEST_NAME ${TESTS_NAMES})
    define_codec_test(${TEST_NAME})
endforeach()
