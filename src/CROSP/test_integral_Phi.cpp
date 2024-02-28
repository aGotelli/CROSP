#include "CROSP/CROSP/cosserat_rod.hpp"


Eigen::MatrixXd getIntegralPhi(const double t_X)
{
    const double end = 1;
    const double begin = 0;
    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    const double y = ( 2 * t_X - ( end + begin) ) / ( end - begin );

    constexpr unsigned int n = 5;


    Eigen::VectorXd integral_phi(n);
    integral_phi << y + 1,
                    0.5 * y*y - 0.5,
                    0.5*(y*y*y - y),
                    0.25*(5*y*y*y*y/2.0 - 3*y*y) + 0.125,
                    0.125*(7*y*y*y*y*y -10*y*y*y + 3*y);


    Eigen::MatrixXd integral_Phi = Eigen::MatrixXd::Zero(3, 3*n);
    for(unsigned int i=0; i<3; i++)
        integral_Phi.block<1, n>(i, i*n) = 0.5*integral_phi;
    /*


                WHY WE NEED TO MULTIPLY FOR 0.5??????


     */

    return integral_Phi;
}


Eigen::MatrixXd getPhi(const double t_X)
{
    const double end = 1;
    const double begin = 0;
    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    const double y = ( 2 * t_X - ( end + begin) ) / ( end - begin );

    constexpr unsigned int n = 5;

    Eigen::VectorXd phi(n);
    phi << 1,
            y,
            0.5*(3*y*y - 1),
            0.5*(5*y*y*y - 3*y),
            0.125*(35*y*y*y*y - 30*y*y + 3);


    Eigen::MatrixXd integral_Phi = Eigen::MatrixXd::Zero(3, 3*n);
    for(unsigned int i=0; i<3; i++)
        integral_Phi.block<1, n>(i, i*n) = phi;

    return integral_Phi;
}




int main(int argc, char *argv[])
{
    static constexpr std::array<bool, 6> admitted_deformations = {
        true,
        true,
        true,

        false,
        false,
        false
    };

    static constexpr unsigned int na = std::count(admitted_deformations.begin(),
                                                  admitted_deformations.end(),
                                                  true);

    unsigned int ne = 5;


    CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);



    const unsigned int n = polynomial_representation.getCoordinatesDimension();



    typedef boost::numeric::odeint::runge_kutta_dopri5< Eigen::MatrixXd, double,
                                                         Eigen::MatrixXd, double,
                                                         boost::numeric::odeint::vector_space_algebra> Ke_stepper;

    const double X0 = 0.0;
    const double X1 = 0.3;
    const double dX = 0.0005;

    Eigen::MatrixXd Phi_integral = Eigen::MatrixXd::Zero(na, n);

    boost::numeric::odeint::integrate_adaptive(Ke_stepper(), [&](const Eigen::MatrixXd &, Eigen::MatrixXd &t_dKeeds, const double t_X){
        const auto Phi = polynomial_representation.getPhi( t_X );

        t_dKeeds = Phi;
    }, Phi_integral, X0, X1, dX);


    std::cout << "Phi_integral : \n" << Phi_integral << "\n\n";


    Phi_integral.setZero();
    boost::numeric::odeint::integrate_adaptive(Ke_stepper(), [&](const Eigen::MatrixXd &, Eigen::MatrixXd &t_dKeeds, const double t_X){
        const auto Phi = getPhi( t_X );

        t_dKeeds = Phi;
    }, Phi_integral, X0, X1, dX);


    std::cout << "Phi_integral with my fun : \n" << Phi_integral << "\n\n";

    std::cout << "Phi_integral analytic : \n" << getIntegralPhi(X1) << "\n\n";

    return 0;
}
