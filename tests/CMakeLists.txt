




include_directories(../src/include)
include_directories(include)






#add_executable(test_idm test_idm.cpp)
#target_link_libraries(test_idm
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark

#)
#target_compile_definitions(test_idm
#    PRIVATE
#        DEVELOPER
#)



#add_executable(test_tidm test_tidm.cpp)
#target_link_libraries(test_tidm
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark

#)
#target_compile_definitions(test_tidm
#    PRIVATE
#        DEVELOPER
#)




#add_executable(benchmark_forward_kinematics benchmark_forward_kinematics.cpp)
#target_link_libraries(benchmark_forward_kinematics PRIVATE ${PROJECT_NAME} benchmark::benchmark Eigen3::Eigen)








#add_executable(test_different_modes test_different_modes.cpp)
#target_link_libraries(test_different_modes
#    PRIVATE
#        utilities
#        ${PROJECT_NAME}
#)
#target_compile_definitions(test_different_modes
#    PRIVATE
#        DEVELOPER
#)


#add_executable(rod_shape rod_shape.cpp)
#target_link_libraries(rod_shape
#    PRIVATE
#        utilities
#        ${PROJECT_NAME}
#        Qt5::Widgets
#        Qt5::Charts
#)
#target_compile_definitions(rod_shape
#    PRIVATE
#        DEVELOPER
#)

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

add_executable(test_tmp_stepper test_tmp_stepper.cpp)
target_link_libraries(test_tmp_stepper
    PRIVATE
        utilities
        ${PROJECT_NAME}
        benchmark::benchmark
)
target_compile_definitions(test_tmp_stepper
    PRIVATE
        DEVELOPER
)




add_executable(test_with_MATLAB_copy test_with_MATLAB_copy.cpp)
target_link_libraries(test_with_MATLAB_copy
    PRIVATE
        utilities
        ${PROJECT_NAME}
)
target_compile_definitions(test_with_MATLAB_copy
    PRIVATE
        DEVELOPER
)



#add_executable(test_with_MATLAB_runge_kutta test_with_MATLAB_runge_kutta.cpp)
#target_link_libraries(test_with_MATLAB_runge_kutta
#    PRIVATE
#        utilities
#        ${PROJECT_NAME}
#)
#target_compile_definitions(test_with_MATLAB_runge_kutta
#    PRIVATE
#        DEVELOPER
#        TESTING
#)

