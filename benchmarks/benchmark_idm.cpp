#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include <benchmark/benchmark.h>



static constexpr std::array<bool, 6> admitted_deformations = {
    false,
    true,
    false,

    false,
    false,
    false
};

static constexpr unsigned int na = std::count(admitted_deformations.begin(),
                                              admitted_deformations.end(),
                                              true);

constexpr unsigned int number_of_Chebyshev_points = 21;


void benchmarkIDM(::benchmark::State &t_state)
{

    const unsigned int ne = t_state.range(0);

    const unsigned int coordinated_dimension = na * ne;


    t_state.counters = {
      {"na", na},
      {"ne", ne}
    };



    auto polynomial_representation =
            std::make_shared<::CROSP::polynomial_representation::PolynomialRepresentation>(admitted_deformations, ne);

    auto strain_param =
            std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(polynomial_representation, number_of_Chebyshev_points);

    ::CROSP::CosseratRod rod(strain_param);

    ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


    Eigen::VectorXd q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);


    while(t_state.KeepRunning()){

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1);

        const Eigen::Vector3d r = rod.getKinematicsAtTip().m_pose.m_position;

        if(std::isnan(r.x()) || std::isnan(r.y()) || std::isnan(r.z())){
            t_state.SkipWithError("Result is nan!");
        }
    }
};






int main(int argc, char *argv[])
{


    const unsigned int repetitions = 20;


    std::vector<unsigned int> ne_stack = {3, 4, 5, 6};


    const std::string benchmark_name = "IDM";


    for(const auto ne : ne_stack)
        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Arg(ne)->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);







    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
