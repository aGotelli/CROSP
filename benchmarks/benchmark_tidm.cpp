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


void benchmarkTIDM(::benchmark::State &t_state)
{

    const unsigned int ne = t_state.range(0);

    const unsigned int coordinated_dimension = na * ne;



    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation);

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


    while(t_state.KeepRunning()){

        for(unsigned int i=0; i<coordinated_dimension; i++){
            Delta_q.setZero();
            Delta_q[i] = 1;

            Delta_dot_q = 400*Delta_q;
            Delta_ddot_q = 16000*Delta_q;

            rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

            rod.forwardTangentKinematics();
            rod.backwardTangentDynamics(F1.block<3,1>(0, 0),
                                        F1.block<3,1>(3, 0));
        }



    }
};






int main(int argc, char *argv[])
{


    const unsigned int repetitions = 5;


    std::vector<unsigned int> ne_stack = {3, 4, 5, 6};


    const std::string benchmark_name = "IDM_na" + std::to_string(na) + "_ne";


    for(const auto ne : ne_stack)
        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkTIDM)->Arg(ne)->Repetitions(repetitions);







    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
