#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"


#include <benchmark/benchmark.h>

#include "BS_thread_pool.hpp"

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



const unsigned int number_of_Chebyshev_points = 17;

const unsigned int ne = 6;

const unsigned int coordinates_dimension = na * ne;



#include <chrono>
#include <thread>


int main(int argc, char *argv[])
{

    const unsigned int repetitions = 20;





    auto polynomial_representation =
                std::make_shared<::CROSP::polynomial_representation::PolynomialRepresentation>(admitted_deformations, ne);

    auto strain_parameterisation =
            std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(polynomial_representation, number_of_Chebyshev_points);

    auto rod_properties =
            std::make_shared<::CROSP::rod_properties::RodProperties>(strain_parameterisation);

    //  The set of integrators needed for the IDM
    std::shared_ptr<::CROSP::idm_integrators::IDMIntegrators> idm_integrators =
        std::make_shared<::CROSP::idm_integrators::IDMIntegrators>(strain_parameterisation,
                                                                   rod_properties );


    Eigen::VectorXd qe = Eigen::VectorXd::Zero( coordinates_dimension );


    Eigen::VectorXd dot_qe = Eigen::VectorXd::Zero( coordinates_dimension );
    Eigen::VectorXd ddot_qe = Eigen::VectorXd::Zero( coordinates_dimension );


    const Eigen::Vector3d N1 = Eigen::Vector3d::Zero();
    const Eigen::Vector3d C1 = Eigen::Vector3d::Zero();

    strain_parameterisation->updateStacks(qe, dot_qe, ddot_qe);

    idm_integrators->m_quaternion->integrate( Eigen::Vector4d(1, 0, 0, 0) );
    idm_integrators->m_position->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_velocity->integrate( Eigen::Vector3d::Zero() );
    idm_integrators->m_linear_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_acceleration->integrate( Eigen::Vector3d::Zero() );
    idm_integrators->m_linear_acceleration->integrate( Eigen::Vector3d::Zero() );

    //  Map force and couple into local coordinates
    const auto q1 = idm_integrators->m_quaternion->getStateAtPoint(0);
    const auto R1 = Eigen::Quaterniond(q1[0], q1[1], q1[2], q1[3]).toRotationMatrix();
    Eigen::Vector3d force_at_tip_local_coord = R1.transpose()*N1;
    Eigen::Vector3d couple_at_tip_local_coord = R1.transpose()*C1;

    idm_integrators->m_internal_forces->integrate( force_at_tip_local_coord );
    idm_integrators->m_internal_couples->integrate( couple_at_tip_local_coord );

    idm_integrators->m_generalised_forces->integrate( Eigen::VectorXd::Zero(coordinates_dimension) );






    std::vector<std::shared_ptr<::CROSP::tidm_integrators::TIDMIntegrators>> tidm_integrators;

    //tidm_integrators.resize(coordinates_dimension);

    const double a = 400;
    const double b = 160000;



    for(unsigned int i=0; i<coordinates_dimension; i++){
        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinates_dimension);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinates_dimension);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinates_dimension);

        Delta_q[i] = 1;
        Delta_dot_q = a * Delta_q;
        Delta_ddot_q = b * Delta_q;

        std::shared_ptr<::CROSP::strain_parameterisation::StrainParameterisation> strain_parameterisation_Delta =
                std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(strain_parameterisation,
                                                                                           ::LieAlgebra::Vector6d::Zero());

        strain_parameterisation_Delta->updateStacks(Delta_q, Delta_dot_q, Delta_ddot_q);


        tidm_integrators.push_back(
                    std::make_shared<::CROSP::tidm_integrators::TIDMIntegrators>(strain_parameterisation,
                                                                                 strain_parameterisation_Delta,
                                                                                 idm_integrators,
                                                                                 rod_properties)
                    );

    }



    Eigen::MatrixXd Jacobian(coordinates_dimension, coordinates_dimension);




    BS::thread_pool pool( std::thread::hardware_concurrency() );
    std::string name = "Jacobian tidm multithreads na=" + std::to_string(na)+" ne="+std::to_string(ne)+" threads=" +std::to_string(pool.get_thread_count());
    ::benchmark::RegisterBenchmark(name.c_str(), [&](::benchmark::State &t_state){

        //  Map force and couple into local coordinates
        Eigen::Vector3d Delta_force_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_force_at_tip;
        Eigen::Vector3d Delta_couple_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_couple_at_tip;

        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinates_dimension);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinates_dimension);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinates_dimension);

        while(t_state.KeepRunning()){
            for(unsigned int i=0; i<coordinates_dimension; i++){
                pool.push_task([&, i](){


                    //  Integrate Delta zeta
                    tidm_integrators[i]->m_Delta_rotation->solveSystem();
                    tidm_integrators[i]->m_Delta_position->solveSystem();

                    //  Integrate Delta eta
                    tidm_integrators[i]->m_Delta_angular_velocity->solveSystem();
                    tidm_integrators[i]->m_Delta_linear_velocity->solveSystem();

                    //  Integrate Delta dot eta
                    tidm_integrators[i]->m_Delta_angular_acceleration->solveSystem();
                    tidm_integrators[i]->m_Delta_linear_acceleration->solveSystem();


                    tidm_integrators[i]->m_Delta_internal_forces->integrate(Delta_force_at_tip_local_coord);
                    tidm_integrators[i]->m_Delta_internal_couples->integrate(Delta_couple_at_tip_local_coord);

                    tidm_integrators[i]->m_Delta_generalised_forces->integrate( Eigen::VectorXd::Zero(coordinates_dimension) );


                    Delta_q[i] = 1;
                    Delta_dot_q = a * Delta_q;
                    Delta_ddot_q = b * Delta_q;


                    Jacobian.col(i) = rod_properties->m_Kee*Delta_q
                                        + rod_properties->m_Dee*Delta_dot_q
                                        - tidm_integrators[i]->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

                });
            }

            pool.wait_for_tasks();
        }
    })/*->Arg(3)->Arg(4)->Arg(5)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)->Arg(20)*/->Repetitions(20)->Unit(::benchmark::kMicrosecond)->UseRealTime();


    name = "Jacobian tidm multithreads loop na=" + std::to_string(na)+" ne="+std::to_string(ne)+" threads=" +std::to_string(pool.get_thread_count());
    ::benchmark::RegisterBenchmark(name.c_str(), [&](::benchmark::State &t_state){

        //  Map force and couple into local coordinates
        Eigen::Vector3d Delta_force_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_force_at_tip;
        Eigen::Vector3d Delta_couple_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_couple_at_tip;

        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinates_dimension);
        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinates_dimension);
        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinates_dimension);

        while(t_state.KeepRunning()){

            pool.push_loop(coordinates_dimension,
                           [&](const int begin, const int end){

                for(int i=begin; i<end; i++){
                    //  Integrate Delta zeta
                    tidm_integrators[i]->m_Delta_rotation->solveSystem();
                    tidm_integrators[i]->m_Delta_position->solveSystem();

                    //  Integrate Delta eta
                    tidm_integrators[i]->m_Delta_angular_velocity->solveSystem();
                    tidm_integrators[i]->m_Delta_linear_velocity->solveSystem();

                    //  Integrate Delta dot eta
                    tidm_integrators[i]->m_Delta_angular_acceleration->solveSystem();
                    tidm_integrators[i]->m_Delta_linear_acceleration->solveSystem();


                    tidm_integrators[i]->m_Delta_internal_forces->integrate(Delta_force_at_tip_local_coord);
                    tidm_integrators[i]->m_Delta_internal_couples->integrate(Delta_couple_at_tip_local_coord);

                    tidm_integrators[i]->m_Delta_generalised_forces->integrate( Eigen::VectorXd::Zero(coordinates_dimension) );


                    Delta_q[i] = 1;
                    Delta_dot_q = a * Delta_q;
                    Delta_ddot_q = b * Delta_q;


                    Jacobian.col(i) = rod_properties->m_Kee*Delta_q
                                        + rod_properties->m_Dee*Delta_dot_q
                                        - tidm_integrators[i]->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

                }
            });


            pool.wait_for_tasks();
        }
    })->Repetitions(20)->Unit(::benchmark::kMicrosecond)->UseRealTime();




//    std::vector<unsigned int> ne_stack = {2, 4, 3};



//    for(const auto ne : ne_stack)
//        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkTIDM)->Arg(ne)->Repetitions(repetitions);


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
