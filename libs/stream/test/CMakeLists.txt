#---------------------------------------------------------------------------#
# Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

include(CMTest)

cm_test_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}

        ${CMAKE_WORKSPACE_NAME}::block

        Boost::unit_test_framework)

macro(define_stream_test name)
    cm_test(NAME stream_${name}_test SOURCES ${name}.cpp)

    target_include_directories(stream_${name}_test PRIVATE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"

            ${Boost_INCLUDE_DIRS})

    set_target_properties(stream_${name}_test PROPERTIES CXX_STANDARD 14
            CXX_STANDARD_REQUIRED TRUE)

endmacro()

set(TESTS_NAMES
        "chacha"
)

foreach(TEST_NAME ${TESTS_NAMES})
    define_stream_test(${TEST_NAME})
endforeach()

cm_test(NAME stream_chacha20poly1305_test SOURCES ../../modes/test/aead_chacha20poly1305.cpp)

target_include_directories(stream_chacha20poly1305_test PRIVATE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../modes/include>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../mac/include>")

target_link_libraries(stream_chacha20poly1305_test
        ${CMAKE_WORKSPACE_NAME}::mac)

set_target_properties(stream_chacha20poly1305_test PROPERTIES CXX_STANDARD 23
        CXX_STANDARD_REQUIRED TRUE)
