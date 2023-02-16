




include_directories(../src/include)
include_directories(include)




add_executable(tests tests.cpp)
target_link_libraries(tests
    PRIVATE
        ${PROJECT_NAME}
)
target_compile_definitions(tests PRIVATE DEVELOPER)


add_executable(test_strain_parameterisation test_strain_parameterisation.cpp)
target_link_libraries(test_strain_parameterisation
    PRIVATE
        ${PROJECT_NAME}
)
target_compile_definitions(test_strain_parameterisation
    PRIVATE
        DEVELOPER
)


add_executable(test_rod_dynamics_vincent test_rod_dynamics_vincent.cpp)
target_link_libraries(test_rod_dynamics_vincent
    PRIVATE
        ${PROJECT_NAME}
        utilities
#        StrainParameterisation
#        RodProperties
#        IDMIntegrators
#        TIDMIntegrators
#        PolynomialRepresentation
)
target_compile_definitions(test_rod_dynamics_vincent
    PRIVATE
        DEVELOPER
)


add_executable(test_tidm test_tidm.cpp)
target_link_libraries(test_tidm
    PRIVATE
        ${PROJECT_NAME}
        utilities
        benchmark::benchmark
#        StrainParameterisation
#        RodProperties
#        IDMIntegrators
#        TIDMIntegrators
#        PolynomialRepresentation
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
        ${PROJECT_NAME}
        utilities
)
target_compile_definitions(test_cosserat_rod
    PRIVATE
        DEVELOPER
)





add_executable(test_static_solution test_static_solution.cpp)
target_link_libraries(test_static_solution
    PRIVATE
        ${PROJECT_NAME}
        utilities
        ${CERES_LIBRARIES}
        Qt5::Widgets
        Qt5::Charts
)
target_compile_definitions(test_static_solution
    PRIVATE
        DEVELOPER
)



add_executable(test_perturbed_idm test_perturbed_idm.cpp)
target_link_libraries(test_perturbed_idm
    PRIVATE
        ${PROJECT_NAME}
        utilities
        ${CERES_LIBRARIES}
        Qt5::Widgets
        Qt5::Charts
)
target_compile_definitions(test_perturbed_idm
    PRIVATE
        DEVELOPER
)




find_package(NEMI)

add_executable(test_internal_elastic_equilibrium test_internal_elastic_equilibrium.cpp)
target_link_libraries(test_internal_elastic_equilibrium
    PRIVATE
        ${PROJECT_NAME}
        utilities
        ${CERES_LIBRARIES}
        Qt5::Widgets
        Qt5::Charts
        NEMI
)
target_compile_definitions(test_internal_elastic_equilibrium
    PRIVATE
        DEVELOPER
)


add_executable(test_DDCT test_DDCT.cpp
    chebyshev_reconstructor.hpp)
target_link_libraries(test_DDCT
    PRIVATE
        ${PROJECT_NAME}
        utilities
)
target_compile_definitions(test_DDCT PRIVATE DEVELOPER)



add_executable(test_different_modes test_different_modes.cpp)
target_link_libraries(test_different_modes
    PRIVATE
        ${PROJECT_NAME}
)
target_compile_definitions(test_different_modes
    PRIVATE
        DEVELOPER
)
