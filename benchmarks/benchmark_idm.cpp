#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

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

constexpr unsigned int number_of_Chebyshev_points = 21;


void benchmarkIDM(::benchmark::State &t_state)
{

    const unsigned int ne =/* t_state.range(0)*/15;

    const unsigned int coordinated_dimension = na * ne;


    t_state.counters = {
      {"na", na},
      {"ne", ne},
      {"Nc", number_of_Chebyshev_points}
    };



    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);



    ::CROSP::CosseratRod<::CROSP::numerical_integrators::explicit_methods::ExplicitIntegrator<>> rod(polynomial_representation, number_of_Chebyshev_points);

    ::LieAlgebra::Vector6d Lambda_X1 = ::LieAlgebra::Vector6d::Zero();


    Eigen::VectorXd q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);



    while(t_state.KeepRunning()){

//        rod.updateParameterisation(q, dot_q, ddot_q);

//        rod.m_cosserat_rod_integrators->updateParameterisation(q, dot_q, ddot_q);

//        rod.m_cosserat_rod_integrators->forwardKinematics();

//        rod.m_cosserat_rod_integrators->backwardDynamics(Lambda_X1);
        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();

        rod.backwardDynamics(Lambda_X1);

    }
};



struct PositionIntegrator : public OSNI::ODEb {

    PositionIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                       const double &t_upper_integration_limit=1.0f,
                       const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_strain_parameterisation_stack->m_Gamma_stack
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    Eigen::Quaterniond m_quaternion_at_point;

};




int main(int argc, char *argv[])
{


//    const unsigned int repetitions = 20;


//    std::vector<unsigned int> ne_stack = {3, 4, 5, 6};


//    const std::string benchmark_name = "IDM";


//    for(const auto ne : ne_stack)
//        ::benchmark::RegisterBenchmark(benchmark_name.c_str(), benchmarkIDM)->Arg(ne)->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);


    ::benchmark::RegisterBenchmark("IDM", benchmarkIDM)->Unit(::benchmark::kMicrosecond);







    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
