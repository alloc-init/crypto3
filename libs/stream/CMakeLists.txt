#---------------------------------------------------------------------------#
# Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

include(CMConfig)
include(CMDeploy)
include(CMSetupVersion)

cm_project(stream WORKSPACE_NAME ${CMAKE_WORKSPACE_NAME} LANGUAGES C CXX)

include(TargetConfiguration)
include(PlatformConfiguration)

option(CRYPTO3_STREAM_CHACHA "Build with CHACHA stream encryption support" TRUE)
option(CRYPTO3_STREAM_SALSA20 "Build with SALSA20 stream encryption support" TRUE)

cm_setup_version(VERSION 0.1.0 PREFIX ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME})

add_library(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE)

set_target_properties(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} PROPERTIES
                      EXPORT_NAME ${CURRENT_PROJECT_NAME})

target_include_directories(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                           "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"

                           ${Boost_INCLUDE_DIRS})

target_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
                      ${CMAKE_WORKSPACE_NAME}::block
                      Boost::container)

cm_deploy(TARGETS ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
          INCLUDE include
          NAMESPACE ${CMAKE_WORKSPACE_NAME}::)

include(CMTest)
cm_add_test_subdirectory(test)
