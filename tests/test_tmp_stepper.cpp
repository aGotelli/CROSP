
#include "CROSP/CROSP/cosserat_rod.hpp"

#include "utilities/Eigen/eigen_io.hpp"

#include <benchmark/benchmark.h>




using namespace ::CROSP::numerical_integrators::runge_kutta;

int main(int argc, char *argv[])
{
    ::benchmark::Initialize(&argc, argv);

    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation;
    unsigned int t_number_of_Chebyshev_points =31;
    std::shared_ptr<::CROSP::rod_properties::RodProperties> rod_properties =
            std::make_shared<::CROSP::rod_properties::RodProperties>();


    ExplicitIntegrator integrator(polynomial_representation,
                                                                               t_number_of_Chebyshev_points,
                                                                               rod_properties);

    std::cout << integrator.printIntegratorProperties() << std::endl;


    unsigned int ne = polynomial_representation.getCoordinatesDimension();
    Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);


    integrator.updateParameterisation(q, dot_q, ddot_q);



    integrator.forwardKinematics();

    auto Lambda_X1 = ::LieAlgebra::Vector6d::Zero();

    integrator.backwardDynamics(Lambda_X1);

    integrator.getIDMStatesObservations();



    ::benchmark::RegisterBenchmark("Backward", [&](::benchmark::State &t_state){

        while(t_state.KeepRunning())
            integrator.backwardDynamics(Lambda_X1);
    });

    ::benchmark::RegisterBenchmark("Backward observer", [&](::benchmark::State &t_state){

        while(t_state.KeepRunning())
            integrator.getIDMStatesObservations();
    });

    ::benchmark::RunSpecifiedBenchmarks();






    return 0;
}
