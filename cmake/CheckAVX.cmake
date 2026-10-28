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

set(AVX512IFMA_CODE "
  #include <immintrin.h>
  int main() {
    __m256i a = _mm256_set1_epi64x(1);
    a = _mm256_madd52lo_epu64(a, a, a);
    a = _mm256_madd52hi_epu64(a, a, a);
    __m512i b = _mm512_set1_epi64(1);
    b = _mm512_madd52lo_epu64(b, b, b);
    b = _mm512_madd52hi_epu64(b, b, b);
    return _mm256_extract_epi64(a, 0) + _mm_cvtsi128_si32(_mm512_castsi512_si128(b));
  }
")

set(SMALL_PRIME_FFT_AVX2_CODE "
  #include <immintrin.h>
  #if !defined(__x86_64__) || !(defined(__GNUC__) || defined(__clang__)) || \
      !defined(__AVX2__) || !defined(__FMA__)
  #error Small-prime AVX2/FMA kernels are not enabled
  #endif
  int main() {
    __m256d a = _mm256_set1_pd(1.0);
    a = _mm256_fmadd_pd(a, a, a);
    return _mm256_cvtsd_f64(a) != 2.0;
  }
")

set(SMALL_PRIME_FFT_AVX512_CODE "
  #include <immintrin.h>
  #if !defined(__x86_64__) || !(defined(__GNUC__) || defined(__clang__)) || \
      !defined(__AVX2__) || !defined(__FMA__) || !defined(__AVX512F__) || !defined(__AVX512DQ__)
  #error Small-prime AVX-512 kernels are not enabled
  #endif
  int main() {
    __m512d a = _mm512_set1_pd(1.0);
    a = _mm512_fmadd_pd(a, a, a);
    a = _mm512_and_pd(a, a);
    __m256d high = _mm512_extractf64x4_pd(a, 1);
    return _mm256_cvtsd_f64(high) != 2.0;
  }
")

function(check_avx_lang lang type flags)
    set(found FALSE)
    set(selected_flags "")
    set(required_flags_save "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_QUIET TRUE)
    foreach(flag ${flags})
        if(NOT found)
            set(CMAKE_REQUIRED_FLAGS "${required_flags_save} ${flag}")
            # Cache the probe, not a previous selection. Source/flag/toolchain changes must be
            # checked again when an existing build directory is reconfigured.
            string(CONCAT probe_context
                "${${type}_CODE}|${CMAKE_REQUIRED_FLAGS}|${CMAKE_REQUIRED_DEFINITIONS}|"
                "${CMAKE_REQUIRED_INCLUDES}|${CMAKE_REQUIRED_LIBRARIES}|"
                "${CMAKE_REQUIRED_LINK_OPTIONS}|${CMAKE_REQUIRED_LINK_DIRECTORIES}|"
                "${CMAKE_${lang}_COMPILER}|${CMAKE_${lang}_COMPILER_ID}|${CMAKE_${lang}_COMPILER_VERSION}|"
                "${CMAKE_${lang}_COMPILER_TARGET}|${CMAKE_${lang}_FLAGS}|${CMAKE_${lang}_FLAGS_DEBUG}|"
                "${CMAKE_${lang}_FLAGS_RELEASE}|${CMAKE_${lang}_FLAGS_RELWITHDEBINFO}|"
                "${CMAKE_${lang}_FLAGS_MINSIZEREL}|${CMAKE_BUILD_TYPE}|${CMAKE_TRY_COMPILE_CONFIGURATION}|"
                "${CMAKE_TRY_COMPILE_TARGET_TYPE}|${CMAKE_SYSROOT}|${CMAKE_OSX_ARCHITECTURES}|"
                "${CMAKE_OSX_SYSROOT}|${CMAKE_OSX_DEPLOYMENT_TARGET}")
            string(SHA256 probe_key "${probe_context}")
            set(probe_var "${lang}_HAS_${type}_${probe_key}")
            get_property(probe_checked GLOBAL PROPERTY "${probe_var}_CHECKED")
            if(NOT probe_checked)
                # Recheck once per configuration, including after linker/toolchain
                # settings outside the probe's immediate context have changed.
                unset(${probe_var} CACHE)
            endif()
            if("${lang}" STREQUAL "CXX")
                check_cxx_source_compiles("${${type}_CODE}" ${probe_var})
            else()
                check_c_source_compiles("${${type}_CODE}" ${probe_var})
            endif()
            set_property(GLOBAL PROPERTY "${probe_var}_CHECKED" TRUE)
            if(${probe_var})
                set(found TRUE)
                set(selected_flags "${flag}")
            endif()
        endif()
    endforeach()
    set(${lang}_${type}_FOUND ${found} CACHE BOOL "${lang} ${type} support" FORCE)
    set(${lang}_${type}_FLAGS "${selected_flags}" CACHE STRING "${lang} ${type} flags" FORCE)
    set(${lang}_${type}_FOUND ${found} PARENT_SCOPE)
    set(${lang}_${type}_FLAGS "${selected_flags}" PARENT_SCOPE)
    mark_as_advanced(${lang}_${type}_FOUND ${lang}_${type}_FLAGS)
endfunction()

# These bundles match the guards in small_prime_fft_arithmetic.hpp. IFMA is a
# separate algebra capability; neither IFMA nor AVX512F alone establishes them.
macro(check_small_prime_fft)
    check_host_cpu_feature(AVX2 "avx2")
    check_host_cpu_feature(FMA "fma")
    check_host_cpu_feature(AVX512F "avx512f")
    check_host_cpu_feature(AVX512DQ "avx512dq")

    foreach(feature SMALL_PRIME_FFT_AVX2 SMALL_PRIME_FFT_AVX512)
        set(CXX_${feature}_FOUND FALSE)
        set(CXX_${feature}_FLAGS "")
    endforeach()
    if(HOST_CPU_AVX2_FOUND AND HOST_CPU_FMA_FOUND)
        # Always export explicit flags, even if the configuring project already supplies them.
        check_avx_lang(CXX SMALL_PRIME_FFT_AVX2 "-mavx2 -mfma")
        if(CXX_SMALL_PRIME_FFT_AVX2_FOUND AND HOST_CPU_AVX512F_FOUND AND HOST_CPU_AVX512DQ_FOUND)
            check_avx_lang(CXX SMALL_PRIME_FFT_AVX512 "-mavx2 -mfma -mavx512f -mavx512dq")
        endif()
    endif()
    foreach(feature SMALL_PRIME_FFT_AVX2 SMALL_PRIME_FFT_AVX512)
        set(CXX_${feature}_FOUND ${CXX_${feature}_FOUND} CACHE BOOL "CXX ${feature} support" FORCE)
        set(CXX_${feature}_FLAGS "${CXX_${feature}_FLAGS}" CACHE STRING "CXX ${feature} flags" FORCE)
        mark_as_advanced(CXX_${feature}_FOUND CXX_${feature}_FLAGS)
    endforeach()
endmacro()

# Both widths use the same host contract; VL is required by the four-lane kernel.
macro(check_avx512ifma)
    # Compiler AVX-512 options can also enable AVX2/FMA implicitly. Check that
    # baseline without depending on the FFT probes or requiring AVX512DQ.
    check_host_cpu_feature(AVX2 "avx2")
    check_host_cpu_feature(FMA "fma")
    check_host_cpu_feature(AVX512F "avx512f")
    check_host_cpu_feature(AVX512IFMA "avx512ifma")
    check_host_cpu_feature(AVX512VL "avx512vl")
    if(HOST_CPU_AVX2_FOUND AND HOST_CPU_FMA_FOUND AND HOST_CPU_AVX512F_FOUND AND
       HOST_CPU_AVX512IFMA_FOUND AND HOST_CPU_AVX512VL_FOUND)
        check_avx_lang(CXX "AVX512IFMA" "-mavx512f -mavx512vl -mavx512ifma;/arch:AVX512")
    else()
        set(CXX_AVX512IFMA_FOUND FALSE)
        set(CXX_AVX512IFMA_FLAGS "")
        set(CXX_AVX512IFMA_FOUND FALSE CACHE BOOL "CXX AVX512IFMA support" FORCE)
        set(CXX_AVX512IFMA_FLAGS "" CACHE STRING "CXX AVX512IFMA flags" FORCE)
    endif()
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
        set(C_AVX_FLAGS "" CACHE STRING "C AVX flags" FORCE)
        set(CXX_AVX_FLAGS "" CACHE STRING "CXX AVX flags" FORCE)
    endif()

    if(HOST_CPU_AVX2_FOUND)
        check_avx_lang(C "AVX2" " ;-mavx2;/arch:AVX2")
        check_avx_lang(CXX "AVX2" " ;-mavx2;/arch:AVX2")
    else()
        set(C_AVX2_FOUND FALSE CACHE BOOL "C AVX2 support" FORCE)
        set(CXX_AVX2_FOUND FALSE CACHE BOOL "CXX AVX2 support" FORCE)
        set(C_AVX2_FLAGS "" CACHE STRING "C AVX2 flags" FORCE)
        set(CXX_AVX2_FLAGS "" CACHE STRING "CXX AVX2 flags" FORCE)
    endif()

    if(HOST_CPU_AVX512F_FOUND)
        check_avx_lang(C "AVX512" " ;-mavx512f;/arch:AVX512")
        check_avx_lang(CXX "AVX512" " ;-mavx512f;/arch:AVX512")
    else()
        set(C_AVX512_FOUND FALSE CACHE BOOL "C AVX512 support" FORCE)
        set(CXX_AVX512_FOUND FALSE CACHE BOOL "CXX AVX512 support" FORCE)
        set(C_AVX512_FLAGS "" CACHE STRING "C AVX512 flags" FORCE)
        set(CXX_AVX512_FLAGS "" CACHE STRING "CXX AVX512 flags" FORCE)
    endif()
endmacro()
