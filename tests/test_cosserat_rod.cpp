#include "CROSP/CROSP/cosserat_rod.hpp"

#include "utilities/Eigen/eigen_io.hpp"





int main()
{
    static constexpr unsigned int number_of_Chebyshev_points = 17;


    static constexpr std::array<bool, 6> admitted_deformations = {
        true,
        true,
        false,

        false,
        false,
        false
    };

    static constexpr unsigned int na = std::count(admitted_deformations.begin(),
                                                  admitted_deformations.end(),
                                                  true);

    static constexpr unsigned int ne = 4;



    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation);


    rod.forwardKinematics();

    Eigen::Vector3d zeros = Eigen::Vector3d::Zero();
    rod.backwardDynamics(zeros, zeros);


    return 0;
}
