




include_directories(../src/include)
include_directories(include)



add_executable(test_with_MATLAB test_with_MATLAB.cpp)
target_link_libraries(test_with_MATLAB
    PRIVATE
        utilities
        ${PROJECT_NAME}
)
target_compile_definitions(test_with_MATLAB
    PRIVATE
        DEVELOPER
)


add_executable(test_with_MATLAB_runge_kutta test_with_MATLAB_runge_kutta.cpp)
target_link_libraries(test_with_MATLAB_runge_kutta
    PRIVATE
        utilities
        ${PROJECT_NAME}
)
target_compile_definitions(test_with_MATLAB_runge_kutta
    PRIVATE
        DEVELOPER
        TESTING
)

