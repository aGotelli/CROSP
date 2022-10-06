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


void benchmarkTIDM(::benchmark::State &t_state)
{

    const unsigned int ne = 3;
    BS::thread_pool pool(t_state.range(0));

    const unsigned int coordinated_dimension = na * ne;

    const unsigned int number_of_Chebyshev_points = 17;



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



    while(t_state.KeepRunning()){


        for(unsigned int i=0; i<coordinated_dimension; i++){

            pool.push_task(

                [&](){

                Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinated_dimension);
                Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
                Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

                //  Variables related the perturbation of the strain parameterisation
                std::shared_ptr<::CROSP::strain_parameterisation::StrainParameterisation> strain_parameterisation_Delta =
                        std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(polynomial_representation,
                                                                                                   ::LieAlgebra::Vector6d::Zero(),
                                                                                                   number_of_Chebyshev_points);

                //  The set of integrators needed for the TIDM
                std::shared_ptr<::CROSP::tidm_integrators::TIDMIntegrators> m_tidm_integrators =
                        std::make_shared<::CROSP::tidm_integrators::TIDMIntegrators>(number_of_Chebyshev_points,
                                                                                    polynomial_representation,
                                                                                    rod.m_strain_parameterisation,
                                                                                    strain_parameterisation_Delta,
                                                                                    rod.m_idm_integrators,
                                                                                    rod.m_rod_properties);


            }

                        );

            Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(coordinated_dimension);
            Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
            Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

            Delta_q.setZero();
            Delta_q[i] = 1;

            Delta_dot_q = 400*Delta_q;
            Delta_ddot_q = 16000*Delta_q;

            pool.push_task(&::CROSP::CosseratRod::TIDM, rod, Delta_q, Delta_dot_q, Delta_ddot_q, F1.block<3,1>(0, 0), F1.block<3,1>(3, 0));

        }

        pool.wait_for_tasks();

    }

};






int main(int argc, char *argv[])
{


    const unsigned int repetitions = 20;


    std::vector<unsigned int> ne_stack = {2, 4, 3};


    const std::string benchmark_name = "TIDM_na" + std::to_string(na) + "_ne";


    for(const auto ne : ne_stack)
        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkTIDM)->Arg(ne)->Repetitions(repetitions);


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
