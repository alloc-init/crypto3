#---------------------------------------------------------------------------#
# Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
#
# SPDX-License-Identifier: MIT
#---------------------------------------------------------------------------#

macro(define_bench_example name)
    add_executable(bench_${name}_example ${name}.cpp)
    target_link_libraries(bench_${name}_example PRIVATE
                          ${CMAKE_WORKSPACE_NAME}_bench
                          ${CMAKE_WORKSPACE_NAME}_algebra
                          ${CMAKE_WORKSPACE_NAME}_multiprecision
                          $<TARGET_NAME_IF_EXISTS:boost_multiprecision>
                      )

    set_target_properties(bench_${name}_example PROPERTIES CXX_STANDARD 14)
endmacro()

set(EXAMPLES_NAMES
    "bls12_381")

foreach(EXAMPLE_NAME ${EXAMPLES_NAMES})
    define_bench_example(${EXAMPLE_NAME})
endforeach()
