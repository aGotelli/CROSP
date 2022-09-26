




include_directories(../src/include)
include_directories(include)





add_executable(benchmark_rod_dynamics benchmark_rod_dynamics.cpp
)
target_link_libraries(benchmark_rod_dynamics
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        math_tools
        utilities
        benchmark::benchmark
)
target_compile_definitions(benchmark_rod_dynamics PRIVATE DEVELOPER)



