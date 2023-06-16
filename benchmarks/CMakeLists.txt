include_directories(../src/include)
include_directories(include)





add_executable(benchmark_idm benchmark_idm.cpp
)
target_link_libraries(benchmark_idm
    PRIVATE
        ${PROJECT_NAME}
        utilities
        benchmark::benchmark
)
target_compile_definitions(benchmark_idm PRIVATE DEVELOPER)

#add_executable(benchmark_tidm benchmark_tidm.cpp
#)
#target_link_libraries(benchmark_tidm
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark
#)
#target_compile_definitions(benchmark_tidm PRIVATE DEVELOPER)


#add_executable(benchmark_magnus benchmark_magnus.cpp
#)
#target_link_libraries(benchmark_magnus
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark
#        optimized OMNI
#        debug OMNI-debug
#)
#target_compile_definitions(benchmark_magnus PRIVATE DEVELOPER)

#add_executable(benchmark_jacobian_tidm benchmark_jacobian_tidm.cpp
#)
#target_link_libraries(benchmark_jacobian_tidm
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark
#)
#target_compile_definitions(benchmark_jacobian_tidm PRIVATE DEVELOPER)




#add_executable(benchmark_idm_components benchmark_idm_components.cpp
#)
#target_link_libraries(benchmark_idm_components
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark
#)
#target_compile_definitions(benchmark_idm_components PRIVATE DEVELOPER)


#add_executable(benchmark_idm_Nc benchmark_idm_Nc.cpp
#)
#target_link_libraries(benchmark_idm_Nc
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark
#)
#target_compile_definitions(benchmark_idm_Nc PRIVATE DEVELOPER)

#add_executable(benchmark_tidm_Nc benchmark_tidm_Nc.cpp
#)
#target_link_libraries(benchmark_tidm_Nc
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark
#)
#target_compile_definitions(benchmark_tidm_Nc PRIVATE DEVELOPER)


#add_executable(benchmark_jacobian_idm_tidm benchmark_jacobian_idm_tidm.cpp
#)
#target_link_libraries(benchmark_jacobian_idm_tidm
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark
#)
#target_compile_definitions(benchmark_jacobian_idm_tidm PRIVATE DEVELOPER)




#add_executable(benchmark_jacobian_tidm_multi_threads
#    benchmark_jacobian_tidm_multi_threads.cpp
#    BS_thread_pool.hpp
#)
#target_link_libraries(benchmark_jacobian_tidm_multi_threads
#    PRIVATE
#        ${PROJECT_NAME}
#        utilities
#        benchmark::benchmark
#)
#target_compile_definitions(benchmark_jacobian_tidm_multi_threads PRIVATE DEVELOPER)

