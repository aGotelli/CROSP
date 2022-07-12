cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(CROSP LANGUAGES CXX VERSION 1.0)
set(ProjectVersion "v_1_0")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)



find_package(Eigen3 3.4 NO_MODULE REQUIRED)



include_directories(include)


find_package(math_tools REQUIRED)
find_package(OSNI REQUIRED)


