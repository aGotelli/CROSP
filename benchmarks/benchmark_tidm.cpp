#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include <benchmark/benchmark.h>







int main(int argc, char *argv[])
{


    const unsigned int repetitions = 20;

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

    std::string benchmark_name;
    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);












    ne = 4;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);

    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);













    ne = 5;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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


        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);














    ne = 6;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);




























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
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);










    ne = 4;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);

    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);













    ne = 5;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);














    ne = 6;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);



























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
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);










    ne = 4;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);

    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);













    ne = 5;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);














    ne = 6;


    na = std::count(admitted_deformations.begin(),
                    admitted_deformations.end(),
                    true);

    benchmark_name.clear();



    benchmark_name = "IDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

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
        }

    })->Repetitions(repetitions);


    benchmark_name.clear();
    benchmark_name = "TIDM_na" + std::to_string(na) + "_ne" + std::to_string(ne);
    ::benchmark::RegisterBenchmark(benchmark_name.c_str(), [&](::benchmark::State &t_state){

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


        Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
        rod.backwardDynamics(F1.block<3,1>(0, 0),
                             F1.block<3,1>(3, 0));


        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


        while(t_state.KeepRunning()){

            for(unsigned int i=0; i< ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;


                Delta_dot_q = 400 * Delta_q;
                Delta_dot_q = 16000 * Delta_q;

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();
                rod.backwardTangentDynamics(Eigen::Vector3d::Zero(),
                                            Eigen::Vector3d::Zero());

            }

        }

    })->Repetitions(repetitions);


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
