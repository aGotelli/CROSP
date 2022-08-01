cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(CROSP LANGUAGES CXX VERSION 1.0)
set(ProjectVersion "v_1_0")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)



find_package(Eigen3 3.4 NO_MODULE REQUIRED)

find_package(math_tools REQUIRED)

find_package(OSNI REQUIRED)


find_package(benchmark REQUIRED)


include_directories(include)



add_library(${PROJECT_NAME} include/material_properties.hpp
                            include/idm_integrators.hpp
                            include/tidm_integrators.hpp
                            include/strain_parameterisation.hpp
                            #src/strain_parameterisation.cpp
                            include/cosserat_rod.hpp
                            src/cosserat_rod.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC math_tools Eigen3::Eigen OSNI)


add_executable(tests tests/tests.cpp)
target_link_libraries(tests ${PROJECT_NAME} math_tools)
target_compile_definitions(tests PRIVATE DEVELOPER)


add_executable(test_strain_parameterisation tests/test_strain_parameterisation.cpp)
target_link_libraries(test_strain_parameterisation ${PROJECT_NAME} math_tools)
target_compile_definitions(test_strain_parameterisation PRIVATE DEVELOPER)

add_executable(test_tidm tests/test_tidm.cpp)
target_link_libraries(test_tidm PRIVATE ${PROJECT_NAME} math_tools benchmark::benchmark)
target_compile_definitions(test_tidm PRIVATE DEVELOPER)


add_executable(benchmark tests/benchmark.cpp)
target_link_libraries(benchmark PRIVATE ${PROJECT_NAME} benchmark::benchmark Eigen3::Eigen)

