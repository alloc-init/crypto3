#---------------------------------------------------------------------------//
#  MIT License
#
#  Copyright (c) 2025 Mikhail Komarov <nemo@allocin.it>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#---------------------------------------------------------------------------//

include(CMConfig)
include(CMSetupVersion)

include(CMDeploy)
include(CMSetupVersion)

cm_project(pubkey WORKSPACE_NAME ${CMAKE_WORKSPACE_NAME})

cm_find_package(CM)
include(CMDeploy)
include(CMSetupVersion)

list(APPEND ${CURRENT_PROJECT_NAME}_LIBRARIES
        ${CMAKE_WORKSPACE_NAME}::algebra
        ${CMAKE_WORKSPACE_NAME}::random)

option(CRYPTO3_PUBKEY_BLS "Build with BLS public key scheme support" TRUE)
option(CRYPTO3_PUBKEY_ECDSA "Build with ECDSA public key scheme support" TRUE)
option(CRYPTO3_PUBKEY_EDDSA "Build with EDDSA public key scheme support" TRUE)
option(CRYPTO3_PUBKEY_SHAMIR "Build with Shamir secret sharing scheme support" TRUE)
option(CRYPTO3_PUBKEY_FELDMAN "Build with Feldman secret sharing scheme support" TRUE)
option(CRYPTO3_PUBKEY_PEDERSEN "Build with Pedersen secret sharing scheme support" TRUE)
option(CRYPTO3_PUBKEY_WEIGHTED_SHAMIR "Build with weighted Shamir secret sharing scheme support" TRUE)

if(CRYPTO3_PUBKEY_BLS)
    add_definitions(-D${CMAKE_UPPER_WORKSPACE_NAME}_HAS_BLS)
endif()

if(CRYPTO3_PUBKEY_ECDSA)
    add_definitions(-D${CMAKE_UPPER_WORKSPACE_NAME}_HAS_ECDSA)
endif()

if(CRYPTO3_PUBKEY_EDDSA)
    add_definitions(-D${CMAKE_UPPER_WORKSPACE_NAME}_HAS_EDDSA)
endif()

cm_setup_version(VERSION 0.1.0 PREFIX ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME})

add_library(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE)

set_target_properties(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} PROPERTIES
        EXPORT_NAME ${CURRENT_PROJECT_NAME})

target_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
        ${${CURRENT_PROJECT_NAME}_LIBRARIES}
        Boost::container)

target_include_directories(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>

        $<$<BOOL:${Boost_FOUND}>:${Boost_INCLUDE_DIRS}>)

cm_deploy(TARGETS ${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
        INCLUDE include
        NAMESPACE ${CMAKE_WORKSPACE_NAME}::)

cm_add_test_subdirectory(test)

if(BUILD_EXAMPLES)
    add_subdirectory(example)
endif()
