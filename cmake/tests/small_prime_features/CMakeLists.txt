cmake_minimum_required(VERSION 3.22)
project(small_prime_feature_consumer LANGUAGES C CXX)

get_filename_component(CRYPTO3_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
option(FEATURE_CHECK_COMPILE_ONLY "Build an object without running target code" OFF)
set(FEATURE_CHECK_HOST "native" CACHE STRING "native, full, avx2, no-fma, no-dq, no-ifma, or unsupported")

if(NOT FEATURE_CHECK_HOST STREQUAL "native")
    if(NOT FEATURE_CHECK_COMPILE_ONLY OR CMAKE_CROSSCOMPILING)
        message(FATAL_ERROR "Simulated host checks require compile-only mode and a native compiler configuration")
    endif()
    # Isolate host-policy simulation here. The compiler probes and Crypto3 targets
    # below are real; a simulated CPU never authorizes executing the built code.
    list(APPEND CMAKE_MODULE_PATH "${CRYPTO3_SOURCE_DIR}/cmake")
    include(CheckHostCPUFeature)
    set(test_features AVX AVX2 FMA AVX512F AVX512DQ AVX512VL AVX512IFMA)
    if(FEATURE_CHECK_HOST STREQUAL "avx2")
        set(test_features AVX AVX2 FMA)
    elseif(FEATURE_CHECK_HOST STREQUAL "no-fma")
        list(REMOVE_ITEM test_features FMA)
    elseif(FEATURE_CHECK_HOST STREQUAL "no-dq")
        list(REMOVE_ITEM test_features AVX512DQ)
    elseif(FEATURE_CHECK_HOST STREQUAL "no-ifma")
        list(REMOVE_ITEM test_features AVX512IFMA)
    elseif(FEATURE_CHECK_HOST STREQUAL "unsupported")
        set(test_features "")
    elseif(NOT FEATURE_CHECK_HOST STREQUAL "full")
        message(FATAL_ERROR "Unknown simulated host: ${FEATURE_CHECK_HOST}")
    endif()
    function(check_host_cpu_feature feature builtin_name)
        if(feature IN_LIST test_features)
            set(supported TRUE)
        else()
            set(supported FALSE)
        endif()
        set(HOST_CPU_${feature}_FOUND ${supported} CACHE BOOL "Simulated host feature" FORCE)
        set(HOST_CPU_${feature}_FOUND ${supported} PARENT_SCOPE)
    endfunction()
endif()

set(BUILD_TESTS OFF CACHE BOOL "Build only requested checks")
set(BUILD_BENCH_TESTS OFF CACHE BOOL "No profiling runs")
enable_testing()
add_subdirectory("${CRYPTO3_SOURCE_DIR}" crypto3 EXCLUDE_FROM_ALL)

if(FEATURE_CHECK_COMPILE_ONLY)
    add_library(feature_consumer OBJECT feature_consumer.cpp)
else()
    add_executable(feature_consumer feature_consumer.cpp)
    add_test(NAME small_prime_feature_consumer COMMAND feature_consumer)
endif()
target_compile_features(feature_consumer PRIVATE cxx_std_23)
# No instruction flags or feature definitions are supplied by this consumer.
target_link_libraries(feature_consumer PRIVATE crypto3::math)

foreach(feature AVX2 AVX512)
    if(CXX_SMALL_PRIME_FFT_${feature}_FOUND)
        set(expect_${feature} 1)
    else()
        set(expect_${feature} 0)
    endif()
endforeach()
set(expect_IFMA 0)
if(expect_AVX512 AND CXX_AVX512IFMA_FOUND)
    set(expect_IFMA 1)
endif()
target_compile_definitions(feature_consumer PRIVATE
    EXPECT_FFT_AVX2=${expect_AVX2} EXPECT_FFT_AVX512=${expect_AVX512} EXPECT_CRT_IFMA=${expect_IFMA})
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/selection.txt"
    "avx2=${expect_AVX2}\navx512=${expect_AVX512}\ncrt_ifma=${expect_IFMA}\n"
    "avx2_flags=${CXX_SMALL_PRIME_FFT_AVX2_FLAGS}\navx512_flags=${CXX_SMALL_PRIME_FFT_AVX512_FLAGS}\n"
    "ifma_found=${CXX_AVX512IFMA_FOUND}\nifma_flags=${CXX_AVX512IFMA_FLAGS}\n"
    "cross=${CMAKE_CROSSCOMPILING}\ncompiler=${CMAKE_CXX_COMPILER_ID}\n"
    "frontend=${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}\n")
