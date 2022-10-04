#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"


#include <benchmark/benchmark.h>

#include "BS_thread_pool.hpp"

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


void benchmarkTIDM(::benchmark::State &t_state)
{

    const unsigned int ne = t_state.range(0);

    const unsigned int coordinated_dimension = na * ne;

    const unsigned int number_of_Chebyshev_points = 30;


    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(number_of_Chebyshev_points, polynomial_representation);

    ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


    Eigen::VectorXd q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

    rod.updateParameterisation(q, dot_q, ddot_q);

    rod.forwardKinematics();
    rod.backwardDynamics(F1.block<3,1>(0, 0),
                         F1.block<3,1>(3, 0));



    Eigen::MatrixXd J(coordinated_dimension, coordinated_dimension);


    std::vector<::CROSP::tidm_integrators::TIDMIntegrators> tidm_integrators_list;

    std::vector<std::shared_ptr<::CROSP::strain_parameterisation::StrainParameterisation>> Delta_strain_parameterisations_list;


    for(unsigned int i=0; i<coordinated_dimension; i++){
        Delta_strain_parameterisations_list.push_back(
            std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(rod.m_polynomial_representation,
                                                                                       ::LieAlgebra::Vector6d::Zero(),
                                                                                       rod.m_number_of_Chebyshev_points)
        );

    }



    for(unsigned int i=0; i<coordinated_dimension; i++){
        tidm_integrators_list.push_back(::CROSP::tidm_integrators::TIDMIntegrators(rod.m_number_of_Chebyshev_points,
                                                                              rod.m_polynomial_representation,
                                                                              rod.m_strain_parameterisation,
                                                                              Delta_strain_parameterisations_list[i],
                                                                              rod.m_idm_integrators,
                                                                              rod.m_rod_properties));
    }

    ::BS::thread_pool pool(8);

    while(t_state.KeepRunning()){

        for(unsigned int i=0; i<coordinated_dimension; i++){


            pool.push_task([i, &Delta_strain_parameterisations_list, &tidm_integrators_list, coordinated_dimension](){
                Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinated_dimension);
                Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
                Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

                Delta_q.setZero();
                Delta_q[i] = 1;

                Delta_dot_q = 400*Delta_q;
                Delta_ddot_q = 16000*Delta_q;

                Delta_strain_parameterisations_list[i]->updateStacks(Delta_q, Delta_dot_q, Delta_ddot_q);

                //  Integrate Delta zeta
                tidm_integrators_list[i].m_Delta_rotation->solveSystem();
                tidm_integrators_list[i].m_Delta_position->solveSystem();

                //  Integrate Delta eta
                tidm_integrators_list[i].m_Delta_angular_velocity->solveSystem();
                tidm_integrators_list[i].m_Delta_linear_velocity->solveSystem();

                //  Integrate Delta dot eta
                tidm_integrators_list[i].m_Delta_angular_acceleration->solveSystem();
                tidm_integrators_list[i].m_Delta_linear_acceleration->solveSystem();


                tidm_integrators_list[i].m_Delta_internal_forces->integrate(Eigen::Vector3d::Zero());
                tidm_integrators_list[i].m_Delta_internal_couples->integrate(Eigen::Vector3d::Zero());

                tidm_integrators_list[i].m_Delta_generalised_forces->integrate(Eigen::VectorXd::Zero(coordinated_dimension));
            });

        }

        pool.wait_for_tasks();


    }
};






int main(int argc, char *argv[])
{


    const unsigned int repetitions = 20;


    std::vector<unsigned int> ne_stack = {3, 4, 5, 6};


    const std::string benchmark_name = "IDM_na" + std::to_string(na) + "_ne";


    for(const auto ne : ne_stack)
        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkTIDM)->Arg(ne)->Repetitions(repetitions);


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
