#include <iostream>



#include <math_tools/Chebyshev/chebyshev_differentiation.hpp>

#include <utilities/Eigen/eigen_io.hpp>

#include "CROSP/numerical_integrators/hybrid_integrator/idm_integrators/idm_integrators.hpp"
// #include "CROSP/CROSP/cosserat_rod.hpp"

int main(int argc, char *argv[])
{

    //  Rod
    ::CROSP::rod_properties::RodProperties* rod_properties =
        new ::CROSP::rod_properties::RodProperties();

    //  Integration parameters
    const unsigned int Nc = 7;
    const unsigned int state_dim = 3;


    //  Integration points
    const auto forward_points  = ::Chebyshev::defineIntegrationPoints(Nc, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD);
    const auto backward_points = ::Chebyshev::defineIntegrationPoints(Nc, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD);

    //  VSA
    const std::array<bool, 6> admitted_deformations = {
        true,
        true,
        true,
        0,
        0,
        0
    };

    const unsigned int na = std::count(admitted_deformations.begin(), admitted_deformations.end(), true);
    const unsigned int modes = 5;
    const unsigned int ne = modes*na;

    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, modes);

    ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* strain_parameterisation_stack =
        new ::CROSP::strain_parameterisation_stack::StrainParameterisationStack(
        polynomial_representation,
        Nc
        );


    ::CROSP::numerical_integrators::hybrid_integrator::idm_integrators::IDMIntegrators idm(strain_parameterisation_stack,
                                                                                           rod_properties);
    // ::CROSP::numerical_integrators::hybrid_integrator::idm_integrators::IDMIntegrators idm_integrators;


    srand( 7 );
    Eigen::VectorXd q   = Eigen::VectorXd::Random(ne);
    Eigen::VectorXd dq  = Eigen::VectorXd::Random(ne)*10;
    Eigen::VectorXd ddq = Eigen::VectorXd::Random(ne)*100;

    //  Print values
    std::cout << "q: " << q << "\n";
    std::cout << "dq: " << dq << "\n";
    std::cout << "ddq: " << ddq << "\n\n";

    strain_parameterisation_stack->updateStrainParameterisation(q, dq, ddq);

    idm.m_forward_kinematics_integrator.integrateODEs();


    std::cout << idm.m_forward_kinematics_integrator.getODE(0)->getStackAsMatrix() << "\n\n";



    // ::CROSP::CosseratRod<> rod(polynomial_representation,
    //                            Nc,
    //                            rod_properties->m_rod_dimensions,
    //                            rod_properties->m_material_properties);
    // rod.m_rod_properties->m_gravity.setZero();
    // rod.updateParameterisation(q, dq, ddq);
    // rod.forwardKinematics();
    // // rod.backwardDynamics(::LieAlgebra::Vector6d::Zero());

    // const auto Omega_reference = rod.m_cosserat_rod_integrators->m_idm_integrators->m_angular_velocity->getStackAsMatrix();

    // std::cout << Omega_reference << "\n\n";


    delete rod_properties;
    delete strain_parameterisation_stack;

    return 0;
}



// #include <iostream>
// #include <memory> // For std::unique_ptr

// // Base class
// class ODECascade {
// public:
//     virtual ~ODECascade() = default; // Virtual destructor for polymorphism
//     virtual void compute() = 0;      // Pure virtual function
// };

// // Derived class: AngularVelocity
// class AngularVelocity : public ODECascade {
// public:
//     void compute() override {
//         std::cout << "AngularVelocity computation\n";
//     }
// };

// // Another derived class: LinearVelocity (specialization of ODECascade)
// class LinearVelocity : public ODECascade {
// private:
//     ODECascade* dependency; // Raw pointer to another ODECascade (non-owning)
// public:
//     // Constructor accepts a raw pointer
//     LinearVelocity(ODECascade* dep) : dependency(dep) {}

//     void compute() override {
//         std::cout << "LinearVelocity computation using dependency:\n";
//         if (dependency) {
//             dependency->compute(); // Use the provided ODECascade instance
//         } else {
//             std::cout << "No dependency provided\n";
//         }
//     }
// };

// int main() {
//     // Create a unique_ptr to manage AngularVelocity
//     std::unique_ptr<ODECascade> angularVel = std::make_unique<AngularVelocity>();

//     // Create a LinearVelocity instance with a raw pointer to AngularVelocity
//     std::unique_ptr<ODECascade> linearVel = std::make_unique<LinearVelocity>(angularVel.get()); // Non-owning reference

//     // Use objects
//     angularVel->compute();  // Direct computation of AngularVelocity
//     linearVel->compute();    // LinearVelocity computation, using AngularVelocity

//     // angularVel is managed by unique_ptr, so no manual delete is required
//     return 0;
// }
