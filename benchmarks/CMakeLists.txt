




include_directories(../src/include)
include_directories(include)





add_executable(benchmark_idm benchmark_idm.cpp
)
target_link_libraries(benchmark_idm
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        math_tools
        utilities
        benchmark::benchmark
)
target_compile_definitions(benchmark_idm PRIVATE DEVELOPER)




add_executable(benchmark_idm_components benchmark_idm_components.cpp
)
target_link_libraries(benchmark_idm_components
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        math_tools
        utilities
        benchmark::benchmark
)
target_compile_definitions(benchmark_idm_components PRIVATE DEVELOPER)


