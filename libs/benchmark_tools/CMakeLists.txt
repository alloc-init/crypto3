#---------------------------------------------------------------------------#
# Copyright (c) 2024 Vasiliy Olekhov <vasiliy.olekhov@nil.foundation>
#
# SPDX-License-Identifier: MIT
#---------------------------------------------------------------------------#

include(CMConfig)
include(CMSetupVersion)
include(CMDeploy)

cm_project(benchmark_tools WORKSPACE_NAME ${CMAKE_WORKSPACE_NAME} LANGUAGES CXX)
cm_setup_version(VERSION 0.1.0 PREFIX ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME})

add_library(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE)
add_library(${CMAKE_WORKSPACE_NAME}::${CURRENT_PROJECT_NAME} ALIAS ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME})

set_target_properties(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} PROPERTIES
        EXPORT_NAME ${CURRENT_PROJECT_NAME})

target_include_directories(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>)

target_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
        Boost::log
)

include(CMTest)

cm_deploy(TARGETS ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
        INCLUDE include
        NAMESPACE ${CMAKE_WORKSPACE_NAME}::)