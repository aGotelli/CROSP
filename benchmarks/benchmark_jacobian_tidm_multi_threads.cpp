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







#include <chrono>
#include <thread>


int main(int argc, char *argv[])
{




    ::benchmark::RegisterBenchmark("4 loops of wait 1000 us", [&](::benchmark::State &t_state){

        BS::thread_pool pool( t_state.range(0) );
        while(t_state.KeepRunning()){
            for(unsigned int i=0; i<4; i++){
                pool.push_task([](){
                    std::this_thread::sleep_for( std::chrono::microseconds(1000) );

                });
            }

            pool.wait_for_tasks();
        }
    })->Arg(3)->Arg(4)->Arg(5)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)->Arg(20)->Repetitions(20)->Unit(::benchmark::kMicrosecond)->UseRealTime();

    const unsigned int number_of_Chebyshev_points = 17;

    const unsigned int ne = 6;

    const unsigned int coordinates_dimension = na * ne;


    const ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation =
            ::CROSP::polynomial_representation::PolynomialRepresentation(admitted_deformations, ne) ;



    //  The set of rod properties
    std::shared_ptr<::CROSP::rod_properties::RodProperties> rod_properties { std::make_shared<::CROSP::rod_properties::RodProperties>(polynomial_representation) };

    //  All the strain releted variables
    std::shared_ptr<::CROSP::strain_parameterisation::StrainParameterisation> strain_parameterisation { std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(polynomial_representation,
                                                                                                                                                                                   number_of_Chebyshev_points) };

    //  The set of integrators needed for the IDM
    std::shared_ptr<::CROSP::idm_integrators::IDMIntegrators> idm_integrators { std::make_shared<::CROSP::idm_integrators::IDMIntegrators>(number_of_Chebyshev_points,
                                                                                                                                           polynomial_representation,
                                                                                                                                           strain_parameterisation,
                                                                                                                                           rod_properties )};


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
                std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(polynomial_representation,
                                                                                           ::LieAlgebra::Vector6d::Zero(),
                                                                                           number_of_Chebyshev_points);

        strain_parameterisation_Delta->updateStacks(Delta_q, Delta_dot_q, Delta_ddot_q);


        tidm_integrators.push_back(
                    std::make_shared<::CROSP::tidm_integrators::TIDMIntegrators>(number_of_Chebyshev_points,
                                                                                 polynomial_representation,
                                                                                 strain_parameterisation,
                                                                                 strain_parameterisation_Delta,
                                                                                 idm_integrators,
                                                                                 rod_properties)
                    );

    }



    Eigen::MatrixXd Jacobian(coordinates_dimension, coordinates_dimension);


//    ::benchmark::RegisterBenchmark("Jacobian column", [&](::benchmark::State &t_state){

//        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinates_dimension);
//        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinates_dimension);
//        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinates_dimension);


//        //  Map force and couple into local coordinates
//        Eigen::Vector3d Delta_force_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_force_at_tip;
//        Eigen::Vector3d Delta_couple_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_couple_at_tip;


//        while(t_state.KeepRunning()){
//            //  Integrate Delta zeta
//            tidm_integrators[0]->m_Delta_rotation->solveSystem();
//            tidm_integrators[0]->m_Delta_position->solveSystem();

//            //  Integrate Delta eta
//            tidm_integrators[0]->m_Delta_angular_velocity->solveSystem();
//            tidm_integrators[0]->m_Delta_linear_velocity->solveSystem();

//            //  Integrate Delta dot eta
//            tidm_integrators[0]->m_Delta_angular_acceleration->solveSystem();
//            tidm_integrators[0]->m_Delta_linear_acceleration->solveSystem();

//            tidm_integrators[0]->m_Delta_internal_forces->integrate(Delta_force_at_tip_local_coord);
//            tidm_integrators[0]->m_Delta_internal_couples->integrate(Delta_couple_at_tip_local_coord);

//            tidm_integrators[0]->m_Delta_generalised_forces->integrate( Eigen::VectorXd::Zero(coordinates_dimension) );


//            Delta_q[0] = 1;
//            Delta_dot_q = a * Delta_q;
//            Delta_ddot_q = b * Delta_q;


//            Jacobian.col(0) = rod_properties->m_Kee*Delta_q
//                                + rod_properties->m_Dee*Delta_dot_q
//                                - tidm_integrators[0]->m_Delta_generalised_forces->getStateAtPoint(::OSNI::ROD_POSITION::BASE);

//        }

//    });

//    ::benchmark::RegisterBenchmark("Jacobian", [&](::benchmark::State &t_state){

//        Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinates_dimension);
//        Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinates_dimension);
//        Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinates_dimension);

//        //  Map force and couple into local coordinates
//        Eigen::Vector3d Delta_force_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_force_at_tip;
//        Eigen::Vector3d Delta_couple_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_couple_at_tip;


//        while(t_state.KeepRunning()){
//            for(unsigned int i=0; i<coordinates_dimension; i++){

//                //  Integrate Delta zeta
//                tidm_integrators[i]->m_Delta_rotation->solveSystem();
//                tidm_integrators[i]->m_Delta_position->solveSystem();

//                //  Integrate Delta eta
//                tidm_integrators[i]->m_Delta_angular_velocity->solveSystem();
//                tidm_integrators[i]->m_Delta_linear_velocity->solveSystem();

//                //  Integrate Delta dot eta
//                tidm_integrators[i]->m_Delta_angular_acceleration->solveSystem();
//                tidm_integrators[i]->m_Delta_linear_acceleration->solveSystem();



//                tidm_integrators[i]->m_Delta_internal_forces->integrate(Delta_force_at_tip_local_coord);
//                tidm_integrators[i]->m_Delta_internal_couples->integrate(Delta_couple_at_tip_local_coord);

//                tidm_integrators[i]->m_Delta_generalised_forces->integrate( Eigen::VectorXd::Zero(coordinates_dimension) );


//                Delta_q[i] = 1;
//                Delta_dot_q = a * Delta_q;
//                Delta_ddot_q = b * Delta_q;


//                Jacobian.col(i) = rod_properties->m_Kee*Delta_q
//                                    + rod_properties->m_Dee*Delta_dot_q
//                                    - tidm_integrators[i]->m_Delta_generalised_forces->getStateAtPoint(::OSNI::ROD_POSITION::BASE);

//            }
//        }
//    });




//    ::benchmark::RegisterBenchmark("Push tasks", [&](::benchmark::State &t_state){
//        while(t_state.KeepRunning()){
//            for(unsigned int i=0; i<coordinates_dimension; i++){
//                pool.push_task([&, i](){
//                    int a = 0;
//                    a++;
//                });
//            }

//            pool.wait_for_tasks();
//        }
//    });



    std::string name = "Jacobian tidm multithreads na=" + std::to_string(na)+" ne="+std::to_string(ne)+" threads=";
    ::benchmark::RegisterBenchmark(name.c_str(), [&](::benchmark::State &t_state){

        BS::thread_pool pool( t_state.range(0) );
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


                    //  Map force and couple into local coordinates
                    Eigen::Vector3d Delta_force_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_force_at_tip;
                    Eigen::Vector3d Delta_couple_at_tip_local_coord = Eigen::Vector3d::Zero();//tip_pose.getRotationMatrix().transpose()*t_Delta_couple_at_tip;


                    tidm_integrators[i]->m_Delta_internal_forces->integrate(Delta_force_at_tip_local_coord);
                    tidm_integrators[i]->m_Delta_internal_couples->integrate(Delta_couple_at_tip_local_coord);

                    tidm_integrators[i]->m_Delta_generalised_forces->integrate( Eigen::VectorXd::Zero(coordinates_dimension) );


                    Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinates_dimension);
                    Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinates_dimension);
                    Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinates_dimension);

                    Delta_q[i] = 1;
                    Delta_dot_q = a * Delta_q;
                    Delta_ddot_q = b * Delta_q;


                    Jacobian.col(i) = rod_properties->m_Kee*Delta_q
                                        + rod_properties->m_Dee*Delta_dot_q
                                        - tidm_integrators[i]->m_Delta_generalised_forces->getStateAtPoint(::OSNI::ROD_POSITION::BASE);

                });
            }

            pool.wait_for_tasks();
        }
    })->Arg(3)->Arg(4)->Arg(5)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)->Arg(20)->Repetitions(20)->Unit(::benchmark::kMicrosecond)->UseRealTime();

//    const unsigned int repetitions = 20;


//    std::vector<unsigned int> ne_stack = {2, 4, 3};



//    for(const auto ne : ne_stack)
//        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkTIDM)->Arg(ne)->Repetitions(repetitions);


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
