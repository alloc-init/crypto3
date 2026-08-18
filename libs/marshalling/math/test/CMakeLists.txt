include(CMTest)

cm_test_link_libraries(${CMAKE_WORKSPACE_NAME}_${CURRENT_PROJECT_NAME}
                       crypto3::algebra
                       crypto3::math
                       ${CMAKE_WORKSPACE_NAME}::crypto3_algebra
                       ${CMAKE_WORKSPACE_NAME}::core
                       Boost::unit_test_framework)

cm_test(NAME marshalling_math_matrix_test SOURCES matrix.cpp)
set_target_properties(marshalling_math_matrix_test PROPERTIES
                      CXX_STANDARD 23
                      CXX_STANDARD_REQUIRED TRUE)
