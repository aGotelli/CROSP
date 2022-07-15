cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(CROSP LANGUAGES CXX VERSION 1.0)
set(ProjectVersion "v_1_0")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)



find_package(Eigen3 3.4 NO_MODULE REQUIRED)

find_package(math_tools REQUIRED)

find_package(OSNI REQUIRED)


include_directories(include)



add_library(${PROJECT_NAME} include/integrators.hpp
                            include/cosserat_rod.hpp
                            src/cosserat_rod.cpp)
target_link_libraries(${PROJECT_NAME} math_tools Eigen3::Eigen OSNI)


add_executable(tests tests/tests.cpp)
target_link_libraries(tests ${PROJECT_NAME} math_tools)
target_compile_definitions(tests PRIVATE DEVELOPER)
