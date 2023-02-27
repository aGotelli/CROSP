




include_directories(../src/include)
include_directories(include)









add_executable(test_tidm test_tidm.cpp)
target_link_libraries(test_tidm
    PRIVATE
        ${PROJECT_NAME}
        utilities
        benchmark::benchmark

)
target_compile_definitions(test_tidm
    PRIVATE
        DEVELOPER
)




add_executable(benchmark_forward_kinematics benchmark_forward_kinematics.cpp)
target_link_libraries(benchmark_forward_kinematics PRIVATE ${PROJECT_NAME} benchmark::benchmark Eigen3::Eigen)








add_executable(test_different_modes test_different_modes.cpp)
target_link_libraries(test_different_modes
    PRIVATE
        utilities
        ${PROJECT_NAME}
)
target_compile_definitions(test_different_modes
    PRIVATE
        DEVELOPER
)
