#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include <benchmark/benchmark.h>



static constexpr std::array<bool, 6> admitted_deformations = {
    true,
    false,
    false,

    false,
    false,
    false
};

static constexpr unsigned int na = std::count(admitted_deformations.begin(),
                                              admitted_deformations.end(),
                                              true);

static constexpr unsigned int ne = 4;

static constexpr unsigned int coordinated_dimension = na * ne;


void benchmarkIDM(::benchmark::State &t_state)
{

    const unsigned int Nc = t_state.range(0);


    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation, Nc);

    ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


    Eigen::VectorXd q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

    rod.updateParameterisation(q, dot_q, ddot_q);

    rod.forwardKinematics();
    rod.backwardDynamics(F1.block<3,1>(0, 0),
                         F1.block<3,1>(3, 0));

    Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

    Delta_q.setZero();
    Delta_q[0] = 1;

    Delta_dot_q = 400*Delta_q;
    Delta_ddot_q = 16000*Delta_q;


    while(t_state.KeepRunning()){

            rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

            rod.forwardTangentKinematics();
            rod.backwardTangentDynamics(F1.block<3,1>(0, 0),
                                        F1.block<3,1>(3, 0));

    }
};






int main(int argc, char *argv[])
{


    const unsigned int repetitions = 20;


    std::vector<unsigned int> Nc_stack = {10, 15, 20, 25, 30, 35};


    const std::string benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne) + "_Nc";

    for(const auto Nc : Nc_stack)
        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Arg(Nc)->Repetitions(repetitions);







    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
