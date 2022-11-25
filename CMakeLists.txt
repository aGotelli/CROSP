cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(CROSP LANGUAGES CXX VERSION 2.4.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)



find_package(Eigen3 3.4 NO_MODULE REQUIRED)

#find_package(utilities REQUIRED)

find_package(math_tools REQUIRED)

find_package(OSNI REQUIRED)

find_package(Ceres REQUIRED)

find_package(Qt5Charts REQUIRED)

find_package(benchmark QUIET)

find_package(GTest QUIET)

#   External utilities to configure the package
include(GNUInstallDirs)
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})


#   Give different names for debug and Release
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_DEBUG_POSTFIX "-debug")


add_subdirectory(src)

# testing
if(GTest_FOUND)
  message(STATUS "Building tests")
  enable_testing()

  add_subdirectory(tests)

endif(GTest_FOUND)


# benchmark
if(benchmark_FOUND)
  message(STATUS "Building benchmark")
  enable_testing()

  add_subdirectory(benchmarks)

endif(benchmark_FOUND)

#   Build documentation with Doxygen (if possible)
include(cmake/documentation_module.cmake)

#   Install the library using the default routine
include(cmake/installation_module.cmake)
