
#include "CROSP/CROSP/cosserat_rod.hpp"

#include "utilities/Eigen/eigen_io.hpp"




using namespace ::CROSP::numerical_integrators::runge_kutta;

int main(int argc, char *argv[])
{


    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation;
    unsigned int t_number_of_Chebyshev_points =31;
    std::shared_ptr<::CROSP::rod_properties::RodProperties> rod_properties =
            std::make_shared<::CROSP::rod_properties::RodProperties>();


    ExplicitIntegrator<runge_kutta_dopri5> integrator(polynomial_representation,
                                                                               t_number_of_Chebyshev_points,
                                                                               rod_properties);

    std::cout << integrator.printIntegratorProperties() << std::endl;

    return 0;
}
