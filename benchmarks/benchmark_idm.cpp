#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include <benchmark/benchmark.h>



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

//constexpr unsigned int number_of_Chebyshev_points = 21;




void benchmarkIDM(::benchmark::State &t_state)
{


    const unsigned int ne = 3;

    const unsigned int number_of_Chebyshev_points = 31;//t_state.range(0);

    const unsigned int coordinated_dimension = na * ne;


    t_state.counters = {
      {"na", na},
      {"ne", ne},
      {"Nc", number_of_Chebyshev_points}
    };



    ::CROSP::polynomial_representation::PolynomialRepresentation
            polynomial_representation(admitted_deformations,
                                      ne,
                                      ::CROSP::polynomial_representation::legendre_polynomial_base);


    ::CROSP::CosseratRod<> rod(polynomial_representation,
                             number_of_Chebyshev_points);

    //rod.printProperties();

    ::LieAlgebra::Vector6d Lambda_X1 = ::LieAlgebra::Vector6d::Zero();


    Eigen::VectorXd q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);



    while(t_state.KeepRunning()){

        rod.updateParameterisation(q, dot_q, ddot_q);
//        rod.updateParameterisation(q, dot_q, ddot_q);
//        rod.updateParameterisation(q, dot_q, ddot_q);
//        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();

//        rod.backwardDynamics(Lambda_X1);

    }
}





int main(int argc, char *argv[])
{


    const unsigned int repetitions = 5;
    const unsigned int Nc_min = 11;
    const unsigned int Nc_step = 2;
    const unsigned int Nc_max = 35;


    ::benchmark::RegisterBenchmark("IDM", benchmarkIDM)->Unit(::benchmark::kMicrosecond)->Repetitions(repetitions)/*->DenseRange(Nc_min, Nc_max, Nc_step)*/;







    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
