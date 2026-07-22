include_guard(GLOBAL)

include(CheckCXXSourceRuns)

function(check_host_cpu_feature feature builtin_name)
    string(TOUPPER "${feature}" feature_upper)
    set(result_var "HOST_CPU_${feature_upper}_FOUND")

    if(CMAKE_CROSSCOMPILING)
        set(${result_var} FALSE CACHE BOOL "${feature} available on host CPU" FORCE)
        set(${result_var} FALSE PARENT_SCOPE)
        return()
    endif()

    if(feature_upper STREQUAL "SSE2")
        set(msvc_test "max_leaf >= 1 && (leaf1[3] & (1 << 26))")
    elseif(feature_upper STREQUAL "SSE3")
        set(msvc_test "max_leaf >= 1 && (leaf1[2] & (1 << 0))")
    elseif(feature_upper STREQUAL "SSSE3")
        set(msvc_test "max_leaf >= 1 && (leaf1[2] & (1 << 9))")
    elseif(feature_upper STREQUAL "SSE4_1")
        set(msvc_test "max_leaf >= 1 && (leaf1[2] & (1 << 19))")
    elseif(feature_upper STREQUAL "SSE4_2")
        set(msvc_test "max_leaf >= 1 && (leaf1[2] & (1 << 20))")
    elseif(feature_upper STREQUAL "AES")
        set(msvc_test "max_leaf >= 1 && (leaf1[2] & (1 << 25))")
    elseif(feature_upper STREQUAL "AVX")
        set(msvc_test "max_leaf >= 1 && (leaf1[2] & (1 << 27)) && (leaf1[2] & (1 << 28)) && ((_xgetbv(0) & 0x6) == 0x6)")
    elseif(feature_upper STREQUAL "AVX2")
        set(msvc_test "max_leaf >= 7 && (leaf1[2] & (1 << 27)) && (leaf1[2] & (1 << 28)) && ((_xgetbv(0) & 0x6) == 0x6) && (leaf7[1] & (1 << 5))")
    elseif(feature_upper STREQUAL "AVX512F")
        set(msvc_test "max_leaf >= 7 && (leaf1[2] & (1 << 27)) && ((_xgetbv(0) & 0xE6) == 0xE6) && (leaf7[1] & (1 << 16))")
    elseif(feature_upper STREQUAL "BMI2")
        set(msvc_test "max_leaf >= 7 && (leaf7[1] & (1 << 8))")
    elseif(feature_upper STREQUAL "ADX")
        set(msvc_test "max_leaf >= 7 && (leaf7[1] & (1 << 19))")
    else()
        message(FATAL_ERROR "Unsupported host CPU feature: ${feature}")
    endif()

    set(source "
        #if defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
        #include <intrin.h>
        int main() {
            int leaf0[4] = {0, 0, 0, 0};
            int leaf1[4] = {0, 0, 0, 0};
            int leaf7[4] = {0, 0, 0, 0};
            __cpuid(leaf0, 0);
            const int max_leaf = leaf0[0];
            if(max_leaf >= 1) __cpuid(leaf1, 1);
            if(max_leaf >= 7) __cpuidex(leaf7, 7, 0);
            return !(${msvc_test});
        }
        #elif (defined(__GNUC__) || defined(__clang__)) && (defined(__i386__) || defined(__x86_64__))
        int main() {
            __builtin_cpu_init();
            return !__builtin_cpu_supports(\"${builtin_name}\");
        }
        #else
        int main() { return 1; }
        #endif
    ")

    check_cxx_source_runs("${source}" ${result_var})
    set(${result_var} "${${result_var}}" PARENT_SCOPE)
    mark_as_advanced(${result_var})
endfunction()
