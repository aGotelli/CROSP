include_directories(include)

add_library(Base SHARED
    include/${PROJECT_NAME}/base/base.hpp
    base/base.cpp)
target_link_libraries(Base
    PRIVATE
        Eigen3::Eigen
)
target_compile_options(Base
    PRIVATE
        -Werror
        -Wall
        -Wextra
)

add_library(RodProperties SHARED
    include/${PROJECT_NAME}/rod_properties/rod_properties.hpp
    rod_properties/rod_properties.cpp
)
target_link_libraries(RodProperties
    PRIVATE
        Eigen3::Eigen
)
target_compile_options(RodProperties
    PRIVATE
        -Werror
        -Wall
        -Wextra
)

add_library(StrainParameterisation SHARED
        include/${PROJECT_NAME}/strain_parameterisation/strain_parameterisation.hpp
        strain_parameterisation/strain_parameterisation.cpp
)
target_link_libraries(StrainParameterisation
    PRIVATE
        Eigen3::Eigen
)
target_compile_options(RodProperties
    PRIVATE
        -Werror
        -Wall
        -Wextra
)


add_library(IDMIntegrators SHARED
    include/${PROJECT_NAME}/idm_integrators/idm_integrators.hpp
    idm_integrators/idm_integrators.cpp
)
target_link_libraries(IDMIntegrators
    PUBLIC
        OSNI
        math_tools
        Eigen3::Eigen
        StrainParameterisation
        RodProperties
)
target_compile_options(RodProperties
    PRIVATE
        -Werror
        -Wall
        -Wextra
)

add_library(TIDMIntegrators SHARED
    include/${PROJECT_NAME}/tidm_integrators/tidm_integrators.hpp
    tidm_integrators/tidm_integrators.cpp
)
target_link_libraries(TIDMIntegrators
    PUBLIC
        Eigen3::Eigen
        OSNI
        math_tools
        StrainParameterisation
        RodProperties
)
target_compile_options(RodProperties
    PRIVATE
        -Werror
        -Wall
        -Wextra
)



add_library(${PROJECT_NAME} SHARED
                            include/${PROJECT_NAME}/${PROJECT_NAME}/cosserat_rod.hpp
                            ${PROJECT_NAME}/cosserat_rod.cpp
)
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        math_tools
        Eigen3::Eigen
        StrainParameterisation
        RodProperties
        IDMIntegrators
        TIDMIntegrators
        OSNI
)
target_compile_options(RodProperties
    PRIVATE
        -Werror
        -Wall
        -Wextra
)


install(TARGETS ${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

