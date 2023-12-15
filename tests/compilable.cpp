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



    ::CROSP::CosseratRod rod(polynomial_representation,
                             number_of_Chebyshev_points,
                             rod_dimensions.makeCopy(),
                             rod_material);

    rod.printProperties();




    auto rod2 = rod;

    rod2.printProperties();



    Eigen::VectorXd q      =     Eigen::VectorXd::Random(coordinated_dimension);
    Eigen::VectorXd dot_q  =  10*Eigen::VectorXd::Random(coordinated_dimension);
    Eigen::VectorXd ddot_q = 100*Eigen::VectorXd::Random(coordinated_dimension);

    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(coordinated_dimension);

    std::cout << "\n\n\n\t\t\tInitial Value for Lambdas \n\n\n\n";
    std::cout << "rod 1 Lambda(X=0) : \n" << rod.getLambdaAtBase() << "\n";
    std::cout << "rod 2 Lambda(X=0) : \n" << rod2.getLambdaAtBase() << "\n";
    std::cout.flush();


    rod.IDM(q, dot_q, ddot_q);

    std::cout << "\n\n\n\t\t\tIDM rod 1 \n\n\n\n";
    std::cout << "rod 1 Lambda(X=0) : \n" << rod.getLambdaAtBase() << "\n";
    std::cout << "rod 2 Lambda(X=0) : \n" << rod2.getLambdaAtBase() << "\n";
    std::cout.flush();

    rod2.IDM(zeros, zeros, zeros);

    std::cout << "\n\n\n\t\t\tNull IDM rod 2 \n\n\n\n";
    std::cout << "rod 1 Lambda(X=0) : \n" << rod.getLambdaAtBase() << "\n";
    std::cout << "rod 2 Lambda(X=0) : \n" << rod2.getLambdaAtBase() << "\n";
    std::cout.flush();


    rod2.IDM(q, dot_q, ddot_q);

    std::cout << "\n\n\n\t\t\tIDM rod 2 \n\n\n\n";
    std::cout << "rod 1 Lambda(X=0) : \n" << rod.getLambdaAtBase() << "\n";
    std::cout << "rod 2 Lambda(X=0) : \n" << rod2.getLambdaAtBase() << "\n";
    std::cout.flush();

//::CROSP::numerical_integrators::spectral_method::SpectralIntegrators integrator(polynomial_representation,
//                                                                                number_of_Chebyshev_points,
//                                                                                std::make_shared<::CROSP::rod_properties::RodProperties>(rod_dimensions,
//                                                                                                                                rod_material)
//                                                                                  );





//    std::cout << integrator.printIntegratorProperties();
    //auto integrator_copy = integrator;

    //std::is_copy_constructible<::CROSP::numerical_integrators::spectral_method::SpectralIntegrators>();




    std::unique_ptr<::CROSP::CosseratRod<>> rod_uptr =
            std::make_unique<::CROSP::CosseratRod<>>(polynomial_representation,
                                                     number_of_Chebyshev_points,
                                                     rod_dimensions,
                                                     rod_material);




    ContinuumJoint joint(std::move(rod_uptr));



    return 0;
}
