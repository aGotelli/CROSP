cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(CROSP LANGUAGES CXX VERSION 2.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)



find_package(Eigen3 3.4 NO_MODULE REQUIRED)

find_package(utilities REQUIRED)

find_package(math_tools REQUIRED)

find_package(OSNI REQUIRED)

find_package(benchmark REQUIRED)

#   External utilities to configure the package
include(GNUInstallDirs)
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})


#   Give different names for debug and Release
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_DEBUG_POSTFIX "-debug")


include_directories(include)

add_subdirectory(src)


#   Install the library using the default routine
include(cmake/installation_module.cmake)
