# ---------------------------------------------------------------------------#
# Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
# ---------------------------------------------------------------------------#

cm_test_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
    ${CMAKE_WORKSPACE_NAME}::random
)

set(TESTS_NAMES
    "polynomial_dfs_benchmark"
)

foreach(TEST_NAME ${TESTS_NAMES})
    define_math_test(${TEST_NAME})
    set(FULL_NAME actor_math_${TEST_NAME}_test)
    target_compile_options(${FULL_NAME} PRIVATE "-march=native")
    target_link_libraries(${FULL_NAME}
        ${CMAKE_WORKSPACE_NAME}::random
        ${CMAKE_WORKSPACE_NAME}::math
        ${CMAKE_WORKSPACE_NAME}::algebra
        ${CMAKE_WORKSPACE_NAME}::multiprecision
        ${CMAKE_WORKSPACE_NAME}::benchmark_tools
        Boost::unit_test_framework
        Boost::timer
        Boost::random)
endforeach()