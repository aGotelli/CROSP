#include "CROSP/CROSP/cosserat_rod.hpp"


Eigen::MatrixXd getIntegralPhi(const double t_X)
{
    const double end = 1;
    const double begin = 0;
    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    const double y = ( 2 * t_X - ( end + begin) ) / ( end - begin );


    Eigen::VectorXd integral_phi(4);
    integral_phi << y +1,
                    0.5 * y*y - 0.5,
                    0.5*(y*y*y -y),
                    0.5*(5*pow(y, 4)/4.0 - 3*y*y/2.0);


    Eigen::MatrixXd integral_Phi = Eigen::MatrixXd::Zero(3, 3*4);
    for(unsigned int i=0; i<3; i++)
        integral_Phi.block<1, 4>(i, i*4) = integral_phi;

    return integral_Phi;
}




int main(int argc, char *argv[])
{
    static constexpr std::array<bool, 6> admitted_deformations = {
        true,
        true,
        true,

        true,
        false,
        false
    };

    static constexpr unsigned int na = std::count(admitted_deformations.begin(),
                                                  admitted_deformations.end(),
                                                  true);

    const std::vector<unsigned int> ne {
        3, 4, 5, 2
    };

    const unsigned int number_of_Chebyshev_points = 31;

    CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);


    ::CROSP::CosseratRod<::CROSP::numerical_integrators::runge_kutta::RungeKuttaIntegrator> rod(polynomial_representation, number_of_Chebyshev_points);


    std::cout << "Kee : \n" << rod.m_Kee << "\n\n";


    const unsigned int n = rod.m_polynomial_representation.getCoordinatesDimension();

    const Eigen::MatrixXd Ha = rod.m_polynomial_representation.m_B.transpose() * rod.m_rod_properties->m_H * rod.m_polynomial_representation.m_B;


    typedef boost::numeric::odeint::runge_kutta_dopri5< Eigen::MatrixXd, double,
                                                         Eigen::MatrixXd, double,
                                                         boost::numeric::odeint::vector_space_algebra> Ke_stepper;
    Eigen::MatrixXd Kee = Eigen::MatrixXd::Zero(n, n);

    const double X0 = 0.0;
    const double X1 = 1.0;
    const double dX = 0.0005;

    boost::numeric::odeint::integrate_adaptive(Ke_stepper(), [&](const Eigen::MatrixXd &, Eigen::MatrixXd &t_dKeeds, const double t_X){
        const auto Phi = rod.m_polynomial_representation.getPhi( t_X );

        t_dKeeds = Phi.transpose()*Phi;
    }, Kee, X0, X1, dX);

    Kee *= rod.m_rod_properties->m_rod_dimensions.m_L;




    unsigned int index = 0;
    for(unsigned int i=0; i<na; i++){
        unsigned int ne_i = ne[i];
        Kee.block(index, index, ne_i, ne_i) *= Ha(i, i);;

        index += ne_i;
    }

    std::cout << "Kee : \n" << Kee << "\n\n";

    return 0;
}
