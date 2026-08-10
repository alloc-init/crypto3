#---------------------------------------------------------------------------#
# Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
# Copyright (c) 2018-2020 Alexey Moskvin
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#


set(BOOST_CMAKE FALSE)

include(CMDeploy)
include(CMConfig)
include(CMSetupVersion)
include(CheckADX)
include(CheckBMI2)

cm_project(multiprecision WORKSPACE_NAME ${CMAKE_WORKSPACE_NAME} LANGUAGES CXX)

cm_setup_version(VERSION 0.1.0 PREFIX ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME})

add_library(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE)

set_target_properties(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} PROPERTIES
        EXPORT_NAME ${CURRENT_PROJECT_NAME})

target_include_directories(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
        ${Boost_INCLUDE_DIRS})

if(CMAKE_TARGET_ARCHITECTURE STREQUAL "x86_64" OR CMAKE_TARGET_ARCHITECTURE STREQUAL "x86")
        check_bmi2()
        check_adx()

        if(CXX_BMI2_FOUND)
                target_compile_definitions(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
                        "${CMAKE_UPPER_WORKSPACE_NAME}_HAS_BMI2")
                target_compile_options(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE ${CXX_BMI2_FLAGS})
        endif()
        if(CXX_ADX_FOUND)
                target_compile_definitions(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
                        "${CMAKE_UPPER_WORKSPACE_NAME}_HAS_ADX")
                target_compile_options(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE ${CXX_ADX_FLAGS})
        endif()
endif()

set_target_properties(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} PROPERTIES
        EXPORT_NAME ${CURRENT_PROJECT_NAME})

cm_deploy(TARGETS ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
        INCLUDE include
        NAMESPACE ${CMAKE_WORKSPACE_NAME}::)

include(CMTest)
cm_add_test_subdirectory(test)
