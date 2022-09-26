include_directories(include)

add_library(PolynomialRepresentation SHARED
    include/${PROJECT_NAME}/polynomial_representation/polynomial_representation.hpp
    polynomial_representation/polynomial_representation.cpp
)
target_link_libraries(PolynomialRepresentation
    PRIVATE
        Eigen3::Eigen
)
target_compile_options(PolynomialRepresentation
    PRIVATE
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
        PolynomialRepresentation
)
target_compile_options(RodProperties
    PRIVATE
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
        PolynomialRepresentation
)
target_compile_options(RodProperties
    PRIVATE
        -Wall
        -Wextra
)





add_library(${PROJECT_NAME}
    SHARED
        include/${PROJECT_NAME}/${PROJECT_NAME}/cosserat_rod.hpp
        ${PROJECT_NAME}/cosserat_rod.cpp
)
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        Eigen3::Eigen
    PUBLIC
        StrainParameterisation
        RodProperties
        PolynomialRepresentation
)
target_compile_options(RodProperties
    PRIVATE
        -Wall
        -Wextra
)


add_executable(main main.cpp)
target_link_libraries(main ${PROJECT_NAME} Eigen3::Eigen)


install(
    TARGETS
        ${PROJECT_NAME}
    DESTINATION
        ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

