include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CheckHostCPUFeature)

set(AVX_CODE "
  #include <immintrin.h>

  int main()
  {
    __m256 a;
    a = _mm256_set1_ps(0);
    return 0;
  }
")

set(AVX2_CODE "
  #include <immintrin.h>

  int main()
  {
    __m256i a = {0};
    a = _mm256_abs_epi16(a);
    __m256i x;
    _mm256_extract_epi64(x, 0); // we rely on this in our AVX2 code
    return 0;
  }
")

set(AVX512_CODE "
  #include <immintrin.h>

  int main()
  {
    __m512i a = _mm512_set1_epi32(1);
    a = _mm512_add_epi32(a, a);
    return _mm_cvtsi128_si32(_mm512_castsi512_si128(a));
  }
")

macro(check_avx_lang lang type flags)
    set(__FLAG_I 1)
    set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
    foreach(__FLAG ${flags})
        if(NOT ${lang}_${type}_FOUND)
            set(CMAKE_REQUIRED_QUIET TRUE)
            set(CMAKE_REQUIRED_FLAGS ${__FLAG})
            if(lang STREQUAL "CXX")
                check_cxx_source_compiles("${${type}_CODE}" ${lang}_HAS_${type}_${__FLAG_I})
            else()
                check_c_source_compiles("${${type}_CODE}" ${lang}_HAS_${type}_${__FLAG_I})
            endif()
            if(${lang}_HAS_${type}_${__FLAG_I})
                set(${lang}_${type}_FOUND TRUE CACHE BOOL "${lang} ${type} support")
                set(${lang}_${type}_FLAGS "${__FLAG}" CACHE STRING "${lang} ${type} flags")
            endif()
            math(EXPR __FLAG_I "${__FLAG_I}+1")
        endif()
    endforeach()
    set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})
    set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

    if(NOT ${lang}_${type}_FOUND)
        set(${lang}_${type}_FOUND FALSE CACHE BOOL "${lang} ${type} support")
        set(${lang}_${type}_FLAGS "" CACHE STRING "${lang} ${type} flags")
    endif()

    mark_as_advanced(${lang}_${type}_FOUND ${lang}_${type}_FLAGS)

endmacro()

macro(check_avx)
    check_host_cpu_feature(AVX "avx")
    check_host_cpu_feature(AVX2 "avx2")
    check_host_cpu_feature(AVX512F "avx512f")

    if(HOST_CPU_AVX_FOUND)
        check_avx_lang(C "AVX" " ;-mavx;/arch:AVX")
        check_avx_lang(CXX "AVX" " ;-mavx;/arch:AVX")
    else()
        set(C_AVX_FOUND FALSE CACHE BOOL "C AVX support" FORCE)
        set(CXX_AVX_FOUND FALSE CACHE BOOL "CXX AVX support" FORCE)
    endif()

    if(HOST_CPU_AVX2_FOUND)
        check_avx_lang(C "AVX2" " ;-mavx2;/arch:AVX2")
        check_avx_lang(CXX "AVX2" " ;-mavx2;/arch:AVX2")
    else()
        set(C_AVX2_FOUND FALSE CACHE BOOL "C AVX2 support" FORCE)
        set(CXX_AVX2_FOUND FALSE CACHE BOOL "CXX AVX2 support" FORCE)
    endif()

    if(HOST_CPU_AVX512F_FOUND)
        check_avx_lang(C "AVX512" " ;-mavx512f;/arch:AVX512")
        check_avx_lang(CXX "AVX512" " ;-mavx512f;/arch:AVX512")
    else()
        set(C_AVX512_FOUND FALSE CACHE BOOL "C AVX512 support" FORCE)
        set(CXX_AVX512_FOUND FALSE CACHE BOOL "CXX AVX512 support" FORCE)
    endif()
endmacro()
