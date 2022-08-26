




include_directories(../src/include/CROSP)
include_directories(include)




add_executable(tests tests.cpp)
target_link_libraries(tests
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        math_tools
)
target_compile_definitions(tests PRIVATE DEVELOPER)


add_executable(test_strain_parameterisation test_strain_parameterisation.cpp)
target_link_libraries(test_strain_parameterisation
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        math_tools
)
target_compile_definitions(test_strain_parameterisation
    PRIVATE
        DEVELOPER
)


add_executable(test_tidm test_tidm.cpp)
target_link_libraries(test_tidm
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        math_tools
        utilities
        benchmark::benchmark
        StrainParameterisation
        MaterialProperties
        IDMIntegrators
        TIDMIntegrators
)
target_compile_definitions(test_tidm
    PRIVATE
        DEVELOPER
)


add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE ${PROJECT_NAME} benchmark::benchmark Eigen3::Eigen)


add_executable(test_cosserat_rod test_cosserat_rod.cpp)
target_link_libraries(test_cosserat_rod
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        math_tools
        utilities
)
target_compile_definitions(test_cosserat_rod
    PRIVATE
        DEVELOPER
)



