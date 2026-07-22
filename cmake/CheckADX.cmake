include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)

set(ADX_CODE "
  #include <immintrin.h>

  #if !defined(__ADX__) && !defined(_MSC_VER)
  #error ADX is not enabled
  #endif

  int main()
  {
    unsigned int result;
    return _addcarryx_u32(0, 1, 2, &result);
  }
")

macro(check_adx_lang lang flags)
    set(__FLAG_I 1)
    set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})
    foreach(__FLAG ${flags})
        if(NOT ${lang}_ADX_FOUND)
            set(CMAKE_REQUIRED_QUIET TRUE)
            set(CMAKE_REQUIRED_FLAGS ${__FLAG})
            if("${lang}" STREQUAL "CXX")
                check_cxx_source_compiles("${ADX_CODE}" ${lang}_HAS_ADX_${__FLAG_I})
            else()
                check_c_source_compiles("${ADX_CODE}" ${lang}_HAS_ADX_${__FLAG_I})
            endif()
            if(${lang}_HAS_ADX_${__FLAG_I})
                set(${lang}_ADX_FOUND TRUE CACHE BOOL "${lang} ADX support")
                set(${lang}_ADX_FLAGS "${__FLAG}" CACHE STRING "${lang} ADX flags")
            endif()
            math(EXPR __FLAG_I "${__FLAG_I}+1")
        endif()
    endforeach()
    set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})

    if(NOT ${lang}_ADX_FOUND)
        set(${lang}_ADX_FOUND FALSE CACHE BOOL "${lang} ADX support")
        set(${lang}_ADX_FLAGS "" CACHE STRING "${lang} ADX flags")
    endif()

    mark_as_advanced(${lang}_ADX_FOUND ${lang}_ADX_FLAGS)
endmacro()

macro(check_adx)
    check_adx_lang(C " ;-madx")
    check_adx_lang(CXX " ;-madx")
endmacro()
