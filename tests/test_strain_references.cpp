

#include "CROSP/polynomial_representation/polynomial_representation.hpp"
#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"









int main(int argc, char *argv[])
{

    const unsigned int number_of_Chebyshev_points = 6;


    const unsigned int ne = 5;


    const std::array<bool, 6> admitted_deformations = {
        true,
        true,
        true,

        true,
        true,
        false
    };

    const unsigned int na = std::count(admitted_deformations.begin(),
                                       admitted_deformations.end(),
                                       true);

    const unsigned int coordinated_dimension = na * ne;





    ::CROSP::polynomial_representation::PolynomialRepresentation
        polynomial_representation(admitted_deformations,
                                  ne,
                                  //::CROSP::polynomial_representation::chebyshev_polynomial_base
                                  ::CROSP::polynomial_representation::legendre_polynomial_base);

    double radius = 0.01;
    ::CROSP::rod_properties::CircularCrossSectionUPtr circular_cross_section =
        std::make_unique<::CROSP::rod_properties::CircularCrossSection>(radius);

    double width = 0.01;
    double height = 0.02;
    ::CROSP::rod_properties::RectangularCrossSectionUPtr rectangular_cross_section =
        std::make_unique<::CROSP::rod_properties::RectangularCrossSection>(width, height);


    double length = 2.45;
    ::CROSP::rod_properties::RodDimensions rod_dimensions(std::move(rectangular_cross_section),
                                                          length);


    ::CROSP::rod_properties::MaterialProperties rod_material{.m_G=0};
    rod_material.m_E = 741e6;


    ::CROSP::strain_parameterisation_stack::StrainParameterisationStack m_strain_parameterisation_stack(polynomial_representation,
                                                                                              number_of_Chebyshev_points);



    Eigen::VectorXd q = Eigen::VectorXd::Random(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Random(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Random(coordinated_dimension);



    m_strain_parameterisation_stack.updateStrainParameterisation(q, dot_q, ddot_q);


    return 0;
}
