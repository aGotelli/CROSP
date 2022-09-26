#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include <benchmark/benchmark.h>












int main(int argc, char *argv[])
{


    const unsigned int repetitions = 5;

    std::array<bool, 6> admitted_deformations = {
        false,
        true,
        false,

        false,
        false,
        false
    };

    unsigned int ne = 3;


    unsigned int na = std::count(admitted_deformations.begin(),
                                 admitted_deformations.end(),
                                 true);


    auto benchmarkIDM = [&ne, &admitted_deformations](::benchmark::State &t_state)
    {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            rod.updateParameterisation(q, dot_q, ddot_q);

            rod.forwardKinematics();
            rod.backwardDynamics(F1.block<3,1>(0, 0),
                                 F1.block<3,1>(3, 0));

            const Eigen::Vector3d r = rod.getKinematicsAtTip().m_pose.m_position;

            if(std::isnan(r.x()) || std::isnan(r.y()) || std::isnan(r.z())){
                t_state.SkipWithError("Result is nan!");
            }
        }
    };







    std::string benchmark_name;
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 4;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 5;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 6;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);





    admitted_deformations = {
        true,
        true,
        false,

        false,
        false,
        false
    };

    ne = 3;
    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);





    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 4;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 5;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 6;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);







    admitted_deformations = {
        true,
        true,
        true,

        false,
        false,
        false
    };

    ne = 3;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);




    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 4;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 5;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);



    ne = 6;
    benchmark_name.clear();
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Repetitions(repetitions);





    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
