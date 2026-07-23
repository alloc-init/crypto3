#---------------------------------------------------------------------------#
# Copyright (c) 2018-2020 Mikhail Komarov <nemo@nil.foundation>
#
# Distributed under the Boost Software License, Version 1.0
# See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt
#---------------------------------------------------------------------------#

include(CheckHostCPUFeature)
include(CheckCXXSourceCompiles)

set(AES_CODE "
    #include <immintrin.h>
    int main() {
        __m128i value = _mm_setzero_si128();
        value = _mm_aesenc_si128(value, value);
        return _mm_cvtsi128_si32(value);
    }
")

set(SSSE3_CODE "
    #include <tmmintrin.h>
    int main() {
        __m128i value = _mm_setzero_si128();
        value = _mm_abs_epi8(value);
        return _mm_cvtsi128_si32(value);
    }
")

macro(check_sse_compiler_feature feature flags)
    set(flag_index 1)
    set(required_flags_save ${CMAKE_REQUIRED_FLAGS})
    foreach(flag ${flags})
        if(NOT CXX_${feature}_FOUND)
            set(CMAKE_REQUIRED_FLAGS ${flag})
            check_cxx_source_compiles("${${feature}_CODE}" CXX_HAS_${feature}_${flag_index})
            if(CXX_HAS_${feature}_${flag_index})
                set(CXX_${feature}_FOUND TRUE CACHE BOOL "CXX ${feature} support")
                set(CXX_${feature}_FLAGS "${flag}" CACHE STRING "CXX ${feature} flags")
            endif()
            math(EXPR flag_index "${flag_index}+1")
        endif()
    endforeach()
    set(CMAKE_REQUIRED_FLAGS ${required_flags_save})

    if(NOT CXX_${feature}_FOUND)
        set(CXX_${feature}_FOUND FALSE CACHE BOOL "CXX ${feature} support")
        set(CXX_${feature}_FLAGS "" CACHE STRING "CXX ${feature} flags")
    endif()
endmacro()

macro(check_sse)
    check_host_cpu_feature(SSE2 "sse2")
    check_host_cpu_feature(SSE3 "sse3")
    check_host_cpu_feature(SSSE3 "ssse3")
    check_host_cpu_feature(SSE4_1 "sse4.1")
    check_host_cpu_feature(SSE4_2 "sse4.2")
    check_host_cpu_feature(AES "aes")

    if(HOST_CPU_SSSE3_FOUND)
        check_sse_compiler_feature(SSSE3 " ;-mssse3")
    else()
        set(CXX_SSSE3_FOUND FALSE CACHE BOOL "CXX SSSE3 support" FORCE)
        set(CXX_SSSE3_FLAGS "" CACHE STRING "CXX SSSE3 flags" FORCE)
    endif()

    if(HOST_CPU_AES_FOUND)
        check_sse_compiler_feature(AES " ;-maes")
    else()
        set(CXX_AES_FOUND FALSE CACHE BOOL "CXX AES support" FORCE)
        set(CXX_AES_FLAGS "" CACHE STRING "CXX AES flags" FORCE)
    endif()

    set(SSE2_FOUND ${HOST_CPU_SSE2_FOUND} CACHE BOOL "SSE2 available on host CPU" FORCE)
    set(SSE3_FOUND ${HOST_CPU_SSE3_FOUND} CACHE BOOL "SSE3 available on host CPU" FORCE)
    set(SSSE3_FOUND ${CXX_SSSE3_FOUND} CACHE BOOL "SSSE3 available on host CPU and compiler" FORCE)
    set(SSE4_1_FOUND ${HOST_CPU_SSE4_1_FOUND} CACHE BOOL "SSE4.1 available on host CPU" FORCE)
    set(SSE4_2_FOUND ${HOST_CPU_SSE4_2_FOUND} CACHE BOOL "SSE4.2 available on host CPU" FORCE)
    set(AES_FOUND ${CXX_AES_FOUND} CACHE BOOL "AES available on host CPU and compiler" FORCE)

    mark_as_advanced(SSE2_FOUND SSE3_FOUND SSSE3_FOUND SSE4_1_FOUND SSE4_2_FOUND AES_FOUND)
endmacro()
