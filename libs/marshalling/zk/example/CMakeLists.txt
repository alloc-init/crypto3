#---------------------------------------------------------------------------#
# Copyright (c) 2018-2021 Mikhail Komarov <nemo@nil.foundation>
# Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

cm_find_package(CM)
include(CMDeploy)
include(CMSetupVersion)

macro(define_marshalling_example name)
    add_executable(${CMAKE_WORKSPACE_NAME}_${name}_example ${name}.cpp)

    target_link_libraries(${CMAKE_WORKSPACE_NAME}_${name}_example PRIVATE

            crypto3::multiprecision
            crypto3::algebra
            crypto3::zk

            ${CMAKE_WORKSPACE_NAME}::core
            ${CMAKE_WORKSPACE_NAME}::crypto3_multiprecision
            ${CMAKE_WORKSPACE_NAME}::crypto3_algebra
            ${CMAKE_WORKSPACE_NAME}::crypto3_zk

            Boost::container)

    target_include_directories(${CMAKE_WORKSPACE_NAME}_${name}_example PRIVATE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"

            ${Boost_INCLUDE_DIRS})

    set_target_properties(${CMAKE_WORKSPACE_NAME}_${name}_example PROPERTIES
            CXX_STANDARD 20
            CXX_STANDARD_REQUIRED TRUE)
endmacro()

set(EXAMPLES_NAMES
        "r1cs_gg_ppzksnark")

foreach(EXAMPLE_NAME ${EXAMPLES_NAMES})
    define_marshalling_example(${EXAMPLE_NAME})
endforeach()
