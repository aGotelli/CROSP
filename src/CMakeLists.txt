include_directories(include)

#add_library(PolynomialRepresentation SHARED
#    include/${PROJECT_NAME}/polynomial_representation/polynomial_representation.hpp
#    polynomial_representation/polynomial_representation.cpp
#)
#target_link_libraries(PolynomialRepresentation
#    PRIVATE
#        Eigen3::Eigen
#)
#target_compile_options(PolynomialRepresentation
#    PRIVATE
#        -Werror
#        -Wall
#        -Wextra
#)

#add_library(RodProperties SHARED
#    include/${PROJECT_NAME}/rod_properties/rod_properties.hpp
#    rod_properties/rod_properties.cpp
#)
#target_link_libraries(RodProperties
#    PRIVATE
#        Eigen3::Eigen
#        PolynomialRepresentation
#)
#target_compile_options(RodProperties
#    PRIVATE
#        -Werror
#        -Wall
#        -Wextra
#)

#add_library(StrainParameterisation SHARED
#        include/${PROJECT_NAME}/strain_parameterisation/strain_parameterisation.hpp
#        strain_parameterisation/strain_parameterisation.cpp
#)
#target_link_libraries(StrainParameterisation
#    PRIVATE
#        Eigen3::Eigen
#        PolynomialRepresentation
#)
#target_compile_options(RodProperties
#    PRIVATE
#        -Werror
#        -Wall
#        -Wextra
#)


#add_library(IDMIntegrators SHARED
#    include/${PROJECT_NAME}/idm_integrators/idm_integrators.hpp
#    idm_integrators/idm_integrators.cpp
#)
#target_link_libraries(IDMIntegrators
#    PUBLIC
#        OSNI
#        math_tools
#        Eigen3::Eigen
#        StrainParameterisation
#        RodProperties
#        PolynomialRepresentation
#)
#target_compile_options(RodProperties
#    PRIVATE
#        -Werror
#        -Wall
#        -Wextra
#)

#add_library(TIDMIntegrators SHARED
#    include/${PROJECT_NAME}/tidm_integrators/tidm_integrators.hpp
#    tidm_integrators/tidm_integrators.cpp
#)
#target_link_libraries(TIDMIntegrators
#    PUBLIC
#        Eigen3::Eigen
#        OSNI
#        math_tools
#        StrainParameterisation
#        RodProperties
#        PolynomialRepresentation
#)
#target_compile_options(RodProperties
#    PRIVATE
#        -Werror
#        -Wall
#        -Wextra
#)



add_library(${PROJECT_NAME} SHARED
                            include/${PROJECT_NAME}/${PROJECT_NAME}/cosserat_rod.hpp
                            ${PROJECT_NAME}/cosserat_rod.cpp
#                            include/${PROJECT_NAME}/static_initialiser/static_initialiser.hpp
#                            static_initialiser/static_initialiser.cpp


                            include/${PROJECT_NAME}/polynomial_representation/polynomial_representation.hpp
                            polynomial_representation/polynomial_representation.cpp
                            include/${PROJECT_NAME}/rod_properties/rod_properties.hpp
                            rod_properties/rod_properties.cpp
                            include/${PROJECT_NAME}/strain_parameterisation/strain_parameterisation.hpp
                            strain_parameterisation/strain_parameterisation.cpp
#                            include/${PROJECT_NAME}/idm_integrators/idm_integrators.hpp
#                            idm_integrators/idm_integrators.cpp
#                            include/${PROJECT_NAME}/tidm_integrators/tidm_integrators.hpp
#                            tidm_integrators/tidm_integrators.cpp
)
target_link_libraries(${PROJECT_NAME}
    PUBLIC
        math_tools
        Eigen3::Eigen
        OSNI
        ${CERES_LIBRARIES}
#    PUBLIC
#        StrainParameterisation
#        RodProperties
#        IDMIntegrators
#        TIDMIntegrators
#        PolynomialRepresentation
)
target_compile_options(${PROJECT_NAME}
    PRIVATE
        -Wall
        -Wextra
)


add_executable(prove prove.cpp)
target_link_libraries(prove ${PROJECT_NAME})


install(
    TARGETS
        ${PROJECT_NAME}
    DESTINATION
        ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

