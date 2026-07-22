include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CheckHostCPUFeature)

set(BMI2_CODE "
  #include <immintrin.h>

  #if !defined(__BMI2__) && !defined(_MSC_VER)
  #error BMI2 is not enabled
  #endif

  int main()
  {
    return _bzhi_u32(3, 1);
  }
")

macro(check_bmi2_lang lang flags)
    set(__FLAG_I 1)
    set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
    foreach(__FLAG ${flags})
        if(NOT ${lang}_BMI2_FOUND)
            set(CMAKE_REQUIRED_FLAGS ${__FLAG})
            set(CMAKE_REQUIRED_QUIET FALSE)
            if("${lang}" STREQUAL "CXX")
                check_cxx_source_compiles("${BMI2_CODE}" ${lang}_HAS_BMI2_${__FLAG_I})
            else()
                check_c_source_compiles("${BMI2_CODE}" ${lang}_HAS_BMI2_${__FLAG_I})
            endif()
            if(${lang}_HAS_BMI2_${__FLAG_I})
                set(${lang}_BMI2_FOUND TRUE CACHE BOOL "${lang} BMI2 support")
                set(${lang}_BMI2_FLAGS "${__FLAG}" CACHE STRING "${lang} BMI2 flags")
            endif()
            math(EXPR __FLAG_I "${__FLAG_I}+1")
        endif()
    endforeach()
    set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})
    set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

    if(NOT ${lang}_BMI2_FOUND)
        set(${lang}_BMI2_FOUND FALSE CACHE BOOL "${lang} BMI2 support")
        set(${lang}_BMI2_FLAGS "" CACHE STRING "${lang} BMI2 flags")
    endif()

    mark_as_advanced(${lang}_BMI2_FOUND ${lang}_BMI2_FLAGS)
endmacro()

macro(check_bmi2)
    check_host_cpu_feature(BMI2 "bmi2")
    if(HOST_CPU_BMI2_FOUND)
        check_bmi2_lang(C " ;-mbmi2")
        check_bmi2_lang(CXX " ;-mbmi2")
    else()
        set(C_BMI2_FOUND FALSE CACHE BOOL "C BMI2 support" FORCE)
        set(C_BMI2_FLAGS "" CACHE STRING "C BMI2 flags" FORCE)
        set(CXX_BMI2_FOUND FALSE CACHE BOOL "CXX BMI2 support" FORCE)
        set(CXX_BMI2_FLAGS "" CACHE STRING "CXX BMI2 flags" FORCE)
    endif()
endmacro()
