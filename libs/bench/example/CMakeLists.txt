#---------------------------------------------------------------------------#
# Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
#
# SPDX-License-Identifier: MIT
#---------------------------------------------------------------------------#

macro(define_bench_example name)
    set(test_name "bench_${name}_example")
    add_executable(${test_name} ${name}.cpp)
    target_link_libraries(${test_name} PRIVATE
                          ${CMAKE_WORKSPACE_NAME}_bench
                          ${CMAKE_WORKSPACE_NAME}_algebra
                          ${CMAKE_WORKSPACE_NAME}_multiprecision
    Boost::random
                      ${Boost_LIBRARIES}
                      )

    set_target_properties(${test_name} PROPERTIES CXX_STANDARD 17)
endmacro()

set(EXAMPLES_NAMES
    "bls12_381")

foreach(EXAMPLE_NAME ${EXAMPLES_NAMES})
    define_bench_example(${EXAMPLE_NAME})
endforeach()
