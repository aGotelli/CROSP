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

constexpr unsigned int number_of_Chebyshev_points = 17;


void benchmarkTIDM(::benchmark::State &t_state)
{

    const unsigned int ne = /*t_state.range(0)*/5;

    const unsigned int coordinated_dimension = na * ne;



    auto polynomial_representation =
            std::make_shared<::CROSP::polynomial_representation::PolynomialRepresentation>(admitted_deformations, ne);

    auto strain_param =
            std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(polynomial_representation, number_of_Chebyshev_points);

    ::CROSP::CosseratRod rod;

    ::LieAlgebra::Vector6d Lambda_X1 = ::LieAlgebra::Vector6d::Zero();


    Eigen::VectorXd q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

    rod.m_cosserat_rod_integrators->updateParameterisation(q, dot_q, ddot_q);

    rod.m_cosserat_rod_integrators->forwardKinematics();

    rod.m_cosserat_rod_integrators->backwardDynamics(Lambda_X1);



    Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

    ::LieAlgebra::Vector6d Delta_Lambda_X1 = ::LieAlgebra::Vector6d::Zero();

    Delta_q.setZero();
    Delta_q[0] = 1;

    Delta_dot_q = 400*Delta_q;
    Delta_ddot_q = 16000*Delta_q;


    while(t_state.KeepRunning()){

        rod.m_cosserat_rod_integrators->updateDeltaParameterisation(Delta_q, Delta_dot_q, Delta_ddot_q);

        rod.m_cosserat_rod_integrators->forwardTangentKinematics();

        rod.m_cosserat_rod_integrators->backwardTangentDynamics(Delta_Lambda_X1);


    }

    t_state.counters = {
      {"na", na},
      {"ne", ne},
      {"Nc", number_of_Chebyshev_points}
    };
};






int main(int argc, char *argv[])
{


//    const unsigned int repetitions = 20;


//    std::vector<unsigned int> ne_stack = {3, 4, 5, 6};


//    const std::string benchmark_name = "TIDM";


//    for(const auto ne : ne_stack)
//        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkTIDM)->Arg(ne)->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);


    ::benchmark::RegisterBenchmark("TIDM", benchmarkTIDM)->Unit(::benchmark::kMicrosecond)->Repetitions(10);


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
