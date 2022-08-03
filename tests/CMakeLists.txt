









add_executable(tests tests/tests.cpp)
target_link_libraries(tests
    PRIVATE
    ${PROJECT_NAME}
    math_tools
    Eigen3::Eigen)
target_compile_definitions(tests PRIVATE DEVELOPER)


add_executable(test_strain_parameterisation tests/test_strain_parameterisation.cpp)
target_link_libraries(test_strain_parameterisation ${PROJECT_NAME} math_tools)
target_compile_definitions(test_strain_parameterisation PRIVATE DEVELOPER)

add_executable(test_tidm tests/test_tidm.cpp)
target_link_libraries(test_tidm
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        math_tools
        utilities
        benchmark::benchmark
        StrainParameterisation
        MaterialProperties
        IDMIntegrators)
target_compile_definitions(test_tidm PRIVATE DEVELOPER)


add_executable(benchmark tests/benchmark.cpp)
target_link_libraries(benchmark PRIVATE ${PROJECT_NAME} benchmark::benchmark Eigen3::Eigen)



