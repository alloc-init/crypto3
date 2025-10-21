include(CMTest)

cm_find_package(${CMAKE_WORKSPACE_NAME}_random)

cm_test_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
        ${CMAKE_WORKSPACE_NAME}::pkpad
        ${CMAKE_WORKSPACE_NAME}::random

        "$<BOOL:${CRYPTO3_PUBKEY_ELGAMAL}:marshalling::${CMAKE_WORKSPACE_NAME}_zk>"

        Boost::unit_test_framework)

macro(define_pubkey_test name)

    cm_test(NAME pubkey_${name}_test SOURCES ${name}.cpp)

    target_include_directories(pubkey_${name}_test PRIVATE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"

            ${Boost_INCLUDE_DIRS})

    set_target_properties(pubkey_${name}_test PROPERTIES CXX_STANDARD 20)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(pubkey_${name}_test PRIVATE "-fconstexpr-steps=2147483647" "-ftemplate-backtrace-limit=0")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(pubkey_${name}_test PRIVATE "-fconstexpr-ops-limit=4294967295" "-ftemplate-backtrace-limit=0")
    endif()
endmacro()

set(TESTS_NAMES
        "ecdsa"
        "bls"
        "secret_sharing"
        "eddsa"
        "schnorr")

if(CRYPTO3_PUBKEY_ELGAMAL)
#    list(APPEND TESTS_NAMES elgamal_verifiable)
endif()

foreach(TEST_NAME ${TESTS_NAMES})
    define_pubkey_test(${TEST_NAME})
endforeach()
