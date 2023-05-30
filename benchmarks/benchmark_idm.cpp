#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

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

constexpr unsigned int number_of_Chebyshev_points = 17;


void benchmarkIDM(::benchmark::State &t_state)
{

    const unsigned int ne = t_state.range(0);

    const unsigned int coordinated_dimension = na * ne;



    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation);

    ::LieAlgebra::Vector6d F1 = ::LieAlgebra::Vector6d::Zero();


    Eigen::VectorXd q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);


    while(t_state.KeepRunning()){

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();
//        rod.backwardDynamics(F1);

//        const Eigen::Vector3d r = rod.getKinematicsAtTip().m_pose.m_position;

//        if(std::isnan(r.x()) || std::isnan(r.y()) || std::isnan(r.z())){
//            t_state.SkipWithError("Result is nan!");
//        }
    }
};






int main(int argc, char *argv[])
{


//    const unsigned int repetitions = 1;


//    std::vector<unsigned int> ne_stack = {3, 4, 5, 6};


//    const std::string benchmark_name = "IDM_na" + std::to_string(na) + "_ne";


//    for(const auto ne : ne_stack)
//        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Arg(ne)->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);








    ::benchmark::RegisterBenchmark("simple call", [](::benchmark::State &t_state){

        std::array<bool, 6> admitted_deformations = {
                true,
                true,
                true,
                false,
                false,
                false
            };


        std::vector<unsigned int> number_of_modes_stack {
            2,
            3,
            5
        };







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




        std::shared_ptr<::CROSP::numerical_integrators::spectral::SpectralIntegrators> m_integrators {
            std::make_shared<::CROSP::numerical_integrators::spectral::SpectralIntegrators>(strain_parameterisation,
                                                                                   strain_parameterisation_Delta,
                                                                                   polynomial_representation,
                                                                                   m_rod_properties,
                                                                                   number_of_Chebyshev_points)
        };





        const auto ne = polynomial_representation.getCoordinatesDimension();
        Eigen::VectorXd q      = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd dot_q  = Eigen::VectorXd::Zero(ne);
        Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);




        q.setRandom();
        dot_q.setRandom();
        ddot_q.setRandom();



        for(const auto _:t_state){
            m_integrators->updateParameterisation(q, dot_q, ddot_q);
            m_integrators->forwardKinematics();

        }

    })->Unit(::benchmark::kMicrosecond);





    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
