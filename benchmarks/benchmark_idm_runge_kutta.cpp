#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"


#include <benchmark/benchmark.h>



const unsigned int repetitions = 10;




int main(int argc, char *argv[])
{


    const std::string path = "../../../MATLAB/test_different_modes/";

    std::array<bool, 6> admitted_deformations = {
            true,
            true,
            true,
            false,
            false,
            false
        };

    Eigen::VectorXd deformations_stack = Eigen::VectorXd::Zero(6);
    for(unsigned int i=0; i<admitted_deformations.size(); i++)
        deformations_stack[i] = static_cast<unsigned int>( admitted_deformations[i] );

    std::vector<unsigned int> number_of_modes_stack {
        2,
        5,
        5
    };

    const unsigned int number_of_modes = std::accumulate(number_of_modes_stack.begin(),
                                                         number_of_modes_stack.end(),
                                                         0);

    Eigen::VectorXd ne_stack = Eigen::VectorXd::Zero(6);
    unsigned int j =0;
    for(unsigned int i=0; i<admitted_deformations.size(); i++)
        if(admitted_deformations[i])
            ne_stack[i] = number_of_modes_stack[j++];






    const unsigned int number_of_Chebyshev_points = 31;
    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, number_of_modes_stack);


    ::CROSP::strain_parameterisation::StrainParameterisation strain_parameterisation;

    ::CROSP::rod_properties::CircularCrossSection cs;
    const double length = 1.75;
    ::CROSP::rod_properties::RodDimensions rod_dimensions(&cs, length);




    //  The set of rod properties
    ::CROSP::rod_properties::RodProperties m_rod_properties {
        ::CROSP::rod_properties::RodProperties(polynomial_representation, rod_dimensions,
                                               ::CROSP::rod_properties::MaterialProperties())
    };



    //  Variables related the perturbation of the strain parameterisation
    ::CROSP::strain_parameterisation::StrainParameterisation strain_parameterisation_Delta {
        ::CROSP::strain_parameterisation::StrainParameterisation(::LieAlgebra::Vector6d::Zero())
    };




    ::CROSP::numerical_integrators::runge_kutta::RungeKuttaIntegrator m_integrators(strain_parameterisation,
                                                                               strain_parameterisation_Delta,
                                                                               polynomial_representation,
                                                                               m_rod_properties,
                                                                               number_of_Chebyshev_points);

    const auto ne = polynomial_representation.getCoordinatesDimension();
    Eigen::VectorXd q      = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd dot_q  = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);




    ::benchmark::RegisterBenchmark("updateParameterisation",[&](::benchmark::State &t_state){


        t_state.counters = {
            {"number_of_modes", number_of_modes}
        };

        for(const auto _ : t_state)
            m_integrators.updateParameterisation(q, dot_q, ddot_q);
    });


    m_integrators.updateParameterisation(q, dot_q, ddot_q);


    ::benchmark::RegisterBenchmark("forwardKinematics",[&](::benchmark::State &t_state){

        t_state.counters = {
            {"number_of_modes", number_of_modes}
        };

        for(const auto _ : t_state)
            m_integrators.forwardKinematics();
    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);


    m_integrators.forwardKinematics();


    auto zeros = ::LieAlgebra::Vector6d::Zero();

    ::benchmark::RegisterBenchmark("backwardDynamics",[&](::benchmark::State &t_state){

        t_state.counters = {
            {"number_of_modes", number_of_modes}
        };

        for(const auto _ : t_state)
            m_integrators.backwardDynamics(zeros);
    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);

    m_integrators.backwardDynamics(zeros);



    ::benchmark::RegisterBenchmark("forwardTangentKinematics",[&](::benchmark::State &t_state){

        t_state.counters = {
            {"number_of_modes", number_of_modes}
        };

        for(const auto _ : t_state)
            m_integrators.forwardTangentKinematics();

    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);


    ::benchmark::RegisterBenchmark("backwardTangentDynamics",[&](::benchmark::State &t_state){

        t_state.counters = {
            {"number_of_modes", number_of_modes}
        };

        for(const auto _ : t_state)
            m_integrators.backwardTangentDynamics(zeros);

    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
