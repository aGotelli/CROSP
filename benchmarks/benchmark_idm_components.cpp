#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include <benchmark/benchmark.h>












int main(int argc, char *argv[])
{


    const unsigned int repetitions = 5;

    std::array<bool, 6> admitted_deformations = {
        true,
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

    std::string benchmark_name;
    const std::string benchmark_name_base = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);



    benchmark_name = benchmark_name_base + "_quaternion";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);


        while(t_state.KeepRunning()){
            rod.m_idm_integrators->m_quaternion->solveSystem();

            if(std::isnan(rod.m_idm_integrators->m_quaternion->getStateAtPoint(0).x()) ||
                std::isnan(rod.m_idm_integrators->m_quaternion->getStateAtPoint(0).y()) ||
                std::isnan(rod.m_idm_integrators->m_quaternion->getStateAtPoint(0).z()) ||
                std::isnan(rod.m_idm_integrators->m_quaternion->getStateAtPoint(0).w()) )
                t_state.SkipWithError("Result is nan!");
        }


    })->Repetitions(repetitions);



    benchmark_name = benchmark_name_base + "_position";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.m_idm_integrators->m_quaternion->solveSystem();


        while(t_state.KeepRunning())
            rod.m_idm_integrators->m_position->solveSystem();


    })->Repetitions(repetitions);



    benchmark_name = benchmark_name_base + "angular_velocities";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.m_idm_integrators->m_quaternion->solveSystem();
        rod.m_idm_integrators->m_position->solveSystem();

        while(t_state.KeepRunning())
            rod.m_idm_integrators->m_angular_velocity->solveSystem();

    })->Repetitions(repetitions);



    benchmark_name = benchmark_name_base + "linear_velocities";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.m_idm_integrators->m_quaternion->solveSystem();
        rod.m_idm_integrators->m_position->solveSystem();
        rod.m_idm_integrators->m_angular_velocity->solveSystem();

        while(t_state.KeepRunning())
            rod.m_idm_integrators->m_linear_velocity->solveSystem();

    })->Repetitions(repetitions);



    benchmark_name = benchmark_name_base + "angular_accelerations";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.m_idm_integrators->m_quaternion->solveSystem();
        rod.m_idm_integrators->m_position->solveSystem();
        rod.m_idm_integrators->m_angular_velocity->solveSystem();
        rod.m_idm_integrators->m_linear_velocity->solveSystem();

        while(t_state.KeepRunning())
            rod.m_idm_integrators->m_angular_acceleration->solveSystem();


    })->Repetitions(repetitions);



    benchmark_name = benchmark_name_base + "linear_accelerations";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.m_idm_integrators->m_quaternion->solveSystem();
        rod.m_idm_integrators->m_position->solveSystem();
        rod.m_idm_integrators->m_angular_velocity->solveSystem();
        rod.m_idm_integrators->m_linear_velocity->solveSystem();
        rod.m_idm_integrators->m_angular_acceleration->solveSystem();

        while(t_state.KeepRunning())
            rod.m_idm_integrators->m_linear_acceleration->solveSystem();


    })->Repetitions(repetitions);



    benchmark_name = benchmark_name_base + "internal_forces";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.m_idm_integrators->m_quaternion->solveSystem();
        rod.m_idm_integrators->m_position->solveSystem();
        rod.m_idm_integrators->m_angular_velocity->solveSystem();
        rod.m_idm_integrators->m_linear_velocity->solveSystem();
        rod.m_idm_integrators->m_angular_acceleration->solveSystem();

        while(t_state.KeepRunning())
            rod.m_idm_integrators->m_internal_forces->solveSystem();


    })->Repetitions(repetitions);



    benchmark_name = benchmark_name_base + "internal_couples";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.m_idm_integrators->m_quaternion->solveSystem();
        rod.m_idm_integrators->m_position->solveSystem();
        rod.m_idm_integrators->m_angular_velocity->solveSystem();
        rod.m_idm_integrators->m_linear_velocity->solveSystem();
        rod.m_idm_integrators->m_angular_acceleration->solveSystem();

        while(t_state.KeepRunning())
            rod.m_idm_integrators->m_internal_couples->solveSystem();



    })->Repetitions(repetitions);



    benchmark_name = benchmark_name_base + "generalised_forces";
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&ne, &admitted_deformations](::benchmark::State &t_state) {
        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.m_idm_integrators->m_quaternion->solveSystem();
        rod.m_idm_integrators->m_position->solveSystem();
        rod.m_idm_integrators->m_angular_velocity->solveSystem();
        rod.m_idm_integrators->m_linear_velocity->solveSystem();
        rod.m_idm_integrators->m_angular_acceleration->solveSystem();
        rod.m_idm_integrators->m_internal_couples->solveSystem();

        while(t_state.KeepRunning())
            rod.m_idm_integrators->m_generalised_forces->solveSystem();



    })->Repetitions(repetitions);






    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
