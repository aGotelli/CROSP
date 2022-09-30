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


void benchmarkJacobianTIDM(::benchmark::State &t_state)
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

    Eigen::MatrixXd J(coordinated_dimension, coordinated_dimension);

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

            J.col(i) = rod.getTangentStaticEquilibrium(Delta_q);
        }



    }
};


void benchmarkJacobianIDMForward(::benchmark::State &t_state)
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


    Eigen::VectorXd residual = rod.getStaticEquilibrium(q);

    Eigen::VectorXd delta_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd delta_dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd delta_ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

    const double delta = 1e-6;

    Eigen::MatrixXd J(coordinated_dimension, coordinated_dimension);

    while(t_state.KeepRunning()){

        for(unsigned int i=0; i<coordinated_dimension; i++){

            delta_q = q;
            delta_q[i] += delta;

            delta_dot_q = 400*delta_q;
            delta_ddot_q = 16000*delta_q;

            rod.updateParameterisation(delta_q, delta_dot_q, delta_ddot_q);

            rod.forwardKinematics();
            rod.backwardDynamics(F1.block<3,1>(0, 0),
                                 F1.block<3,1>(3, 0));

            const Eigen::VectorXd residual_difference = rod.getStaticEquilibrium(delta_q) - residual;

            J.col(i) = residual_difference / delta;
        }

    }
};



void benchmarkJacobianIDMCentral(::benchmark::State &t_state)
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


    Eigen::VectorXd residual = rod.getStaticEquilibrium(q);

    Eigen::VectorXd delta_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd delta_dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd delta_ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

    const double delta = 1e-6;

    Eigen::MatrixXd J(coordinated_dimension, coordinated_dimension);

    while(t_state.KeepRunning()){

        for(unsigned int i=0; i<coordinated_dimension; i++){

            delta_q = q;
            delta_q[i] += delta;

            delta_dot_q = 400*delta_q;
            delta_ddot_q = 16000*delta_q;

            rod.updateParameterisation(delta_q, delta_dot_q, delta_ddot_q);

            rod.forwardKinematics();
            rod.backwardDynamics(F1.block<3,1>(0, 0),
                                 F1.block<3,1>(3, 0));

            const Eigen::VectorXd increment_residual = rod.getStaticEquilibrium(delta_q);


            delta_q = q;
            delta_q[i] -= delta;

            delta_dot_q = 400*delta_q;
            delta_ddot_q = 16000*delta_q;

            rod.updateParameterisation(delta_q, delta_dot_q, delta_ddot_q);

            rod.forwardKinematics();
            rod.backwardDynamics(F1.block<3,1>(0, 0),
                                 F1.block<3,1>(3, 0));

            const Eigen::VectorXd decrement_residual = rod.getStaticEquilibrium(delta_q);

            const Eigen::VectorXd residual_difference = increment_residual - decrement_residual;

            J.col(i) = residual_difference / delta;
        }

    }
};



int main(int argc, char *argv[])
{


    const unsigned int repetitions = 20;


    std::vector<unsigned int> ne_stack = {3, 4, 5, 6};


    std::string benchmark_name = "Jacobian_TIDM_na" + std::to_string(na) + "_ne";


//    for(const auto ne : ne_stack)
//        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkJacobianTIDM)->Arg(ne)->Repetitions(repetitions);


    benchmark_name = "Jacobian_forward_IDM_na" + std::to_string(na) + "_ne";

    for(const auto ne : ne_stack)
        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkJacobianIDMForward)->Arg(ne)->Repetitions(repetitions);


    benchmark_name = "Jacobian_central_IDM_na" + std::to_string(na) + "_ne";

    for(const auto ne : ne_stack)
        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkJacobianIDMCentral)->Arg(ne)->Repetitions(repetitions);






    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
