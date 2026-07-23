#---------------------------------------------------------------------------#
# Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake"
        "${CMAKE_CURRENT_LIST_DIR}/cmake/packages"
        "${CMAKE_CURRENT_LIST_DIR}/cmake/modules/share/modules/cmake")

set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

include(CMConfig)
include(CMSetupVersion)

cm_workspace(crypto3)
cm_setup_version(VERSION 0.3.0 PREFIX ${CMAKE_WORKSPACE_NAME})

option(BUILD_TESTS "Enable tests" TRUE)# used by CMTest module; If set to false, tests are removed from crypto3::all target. User still can build individual ones.
option(BUILD_BENCH_TESTS "Build performance benchmark tests" FALSE)
option(BUILD_DOCS "Build with configuring Doxygen documentation compiler" FALSE)

# The file compile_commands.json is generated in build directory, so LSP could
# pick it up and guess all include paths, defines and other stuff.
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
            ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-fdiagnostics-color=always)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(-fcolor-diagnostics)
endif()

set(BUILD_WITH_TARGET_ARCHITECTURE "" CACHE STRING "Target build architecture") # TODO: check if we can do cross-compilation
set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_LIST_DIR}/docs" CACHE STRING "Specify doxygen output directory")

include(TargetArchitecture)

if(BUILD_WITH_TARGET_ARCHITECTURE)
    set(CMAKE_TARGET_ARCHITECTURE ${BUILD_WITH_TARGET_ARCHITECTURE})
else()
    target_architecture(CMAKE_TARGET_ARCHITECTURE)
endif()

find_package(Boost REQUIRED COMPONENTS
        container random filesystem log log_setup program_options thread unit_test_framework timer)

function(add_tests)
    if(BUILD_TESTS)
        add_subdirectory(test)
    endif()
endfunction(add_tests)

add_subdirectories("${CMAKE_CURRENT_LIST_DIR}/libs/")
add_subdirectories("${CMAKE_CURRENT_LIST_DIR}/libs/marshalling/")

configure_file(${CMAKE_CURRENT_LIST_DIR}/docs/doxygen/${CMAKE_WORKSPACE_NAME}.doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_WORKSPACE_NAME}.doxyfile @ONLY)