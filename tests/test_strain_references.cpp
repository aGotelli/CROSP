#include "CROSP/CROSP/cosserat_rod.hpp"

#include "utilities/Eigen/eigen_io.hpp"



struct ContinuumJoint {

    ContinuumJoint(std::unique_ptr<::CROSP::CosseratRod<>> t_rod)
        : m_rod(std::move(t_rod))
    {}


    std::unique_ptr<::CROSP::CosseratRod<>> m_rod;
};








int main(int argc, char *argv[])
{

    const unsigned int number_of_Chebyshev_points = 21;


    const unsigned int ne = 5;


    const std::array<bool, 6> admitted_deformations = {
        false,
        true,
        false,

        false,
        false,
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






    return 0;
}
