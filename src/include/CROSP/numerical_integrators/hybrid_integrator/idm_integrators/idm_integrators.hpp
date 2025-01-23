/**
 * \file idm_integrators.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod IDM
 * \date 13-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */



#ifndef IDM_INTEGRATORS_HPP
#define IDM_INTEGRATORS_HPP

#include <memory>

#include "CROSP/rod_properties/rod_properties.hpp"
//#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include "OSNI/OSNI.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"

#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"
#include "ATORS/ATORS.hpp"


/// \brief CROSP::idm_integrators is the namespace containing the integrators for the IDM
namespace CROSP::numerical_integrators::hybrid_integrator::idm_integrators {








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



struct AngularVelocityIntegrator : public OSNI::ODECascade {

    virtual ~AngularVelocityIntegrator(){};

    AngularVelocityIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack *t_strain_parameterisation_stack);

    virtual Eigen::VectorXd computeParameterVectorAtPoint(const unsigned int t_point) final;

    const std::vector<Eigen::Vector3d>* m_dot_K_stack;

};

struct LinearVelocityIntegrator : public OSNI::ODECascade {

    virtual ~LinearVelocityIntegrator() {};

    LinearVelocityIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                             std::shared_ptr<const OSNI::ODECascade> t_angular_velocity_integrator);

    virtual Eigen::VectorXd computeParameterVectorAtPoint(const unsigned int t_point) final;

    const std::vector<Eigen::Vector3d>* m_Gamma_stack;

    const std::vector<Eigen::Vector3d>* m_dot_Gamma_stack;

    std::shared_ptr<const OSNI::ODECascade> m_angular_velocity { nullptr };
};




// struct AngularAccelerationIntegrator : public OSNI::ODEAb {

//     AngularAccelerationIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
//                                   std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
//                                   const double &t_upper_integration_limit=1.0f,
//                                   const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


//     virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


//     virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


//     std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

//     std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
//         m_strain_parameterisation_stack->m_hat_K_stack
//     };

//     std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_K_stack {
//         m_strain_parameterisation_stack->m_hat_dot_K_stack
//     };

//     std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_K_stack {
//         m_strain_parameterisation_stack->m_ddot_K_stack
//     };

//     std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;


// };







// struct LinearAccelerationIntegrator : public OSNI::ODEAb {

//     LinearAccelerationIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
//                                  std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
//                                  std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
//                                  std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
//                                  const double &t_upper_integration_limit=1.0f,
//                                  const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());



//     virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


//     virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


//     std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

//     std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
//         m_strain_parameterisation_stack->m_hat_K_stack
//     };

//     std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_K_stack {
//         m_strain_parameterisation_stack->m_hat_dot_K_stack
//     };

//     std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_Gamma_stack {
//         m_strain_parameterisation_stack->m_hat_Gamma_stack
//     };

//     std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_Gamma_stack{
//         m_strain_parameterisation_stack->m_hat_dot_Gamma_stack
//     };

//     std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_Gamma_stack{
//         m_strain_parameterisation_stack->m_ddot_Gamma_stack
//     };

//     std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;



// };




// struct InternalForcesIntegrator : public OSNI::ODEAb {

//     InternalForcesIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
//                              std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
//                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
//                              std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
//                              std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
//                              std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
//                              std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
//                              const double &t_upper_integration_limit=1.0f,
//                              const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());



//     virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


//     virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


//     virtual Eigen::VectorXd computeDistributedForce(const unsigned int t_point) const;




//     std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

//     std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
//         m_strain_parameterisation_stack->m_hat_K_stack
//     };

//     std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_position;

// };


// struct InternalCouplesIntegrator : public OSNI::ODEAb {

//     InternalCouplesIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
//                               std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
//                               std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
//                               std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
//                               std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
//                               std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
//                               std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
//                               std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
//                               const double &t_upper_integration_limit=1.0f,
//                               const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());




//     virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


//     virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


//     virtual Eigen::VectorXd computeDistributedCouple(const unsigned int) const;


//     std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

//     std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
//         m_strain_parameterisation_stack->m_hat_K_stack
//     };

//     std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
//         m_strain_parameterisation_stack->m_Gamma_stack
//     };

//     std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_position;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces;

// };


// struct GeneralisedForcesIntegrator : public OSNI::ODEb {

//     GeneralisedForcesIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
//                                 std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
//                                 std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
//                                 const double &t_upper_integration_limit=1.0f);

//     virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

//     std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

//     unsigned int m_coordinates_dimension;

//     std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator;
//     std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;

// };







/*!
 * \brief The IDMIntegrators struct contains all the intgrators in order to integrate the kinematics and dynamics of the rod
 *
 * This struct is a collection of integrators to use in order to compute the Forward Kinematics and the Backward Dynamics
 */
struct IDMIntegrators {

    IDMIntegrators(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                   const rod_properties::RodProperties* t_rod_properties);


    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    //  Instance of the rod properties
    const rod_properties::RodProperties* m_rod_properties;


    unsigned int m_states_dimension { 3 };

    std::function<Eigen::MatrixXd(unsigned int)> m_A_function {
        [this](unsigned int t_point) {
            return -::LieAlgebra::skew( m_strain_parameterisation_stack->m_K_stack->at(t_point) );
        }
    };

    OSNI::InCascateODEsIntegrator m_forward_kinematics_integrator {
        OSNI::InCascateODEsIntegrator(::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                                      m_strain_parameterisation_stack->m_number_of_points,
                                      m_states_dimension,
                                      m_A_function)
    };


    // //  Integrator for the quaternions
    // std::unique_ptr<OSNI::ODECascade> m_quaternion {
    //     std::make_shared<QuaternionIntegrator>(m_strain_parameterisation_stack,
    //                                            m_rod_properties->m_rod_dimensions.m_L)
    // };

    // //  Integrator for the positions
    // std::shared_ptr<OSNI::ODESolverInterface> m_position {
    //     std::make_shared<PositionIntegrator>(m_strain_parameterisation_stack,
    //                                          m_quaternion,
    //                                          m_rod_properties->m_rod_dimensions.m_L)
    // };

    // //  Integrator for the angular velocities
    // std::shared_ptr<OSNI::ODESolverInterface> m_angular_velocity {
    //     std::make_shared<AngularVelocityIntegrator>(m_strain_parameterisation_stack,
    //                                                 m_rod_properties->m_rod_dimensions.m_L)
    // };

    // //  Integrator for the linear velocities
    // std::shared_ptr<OSNI::ODESolverInterface> m_linear_velocity {
    //     std::make_shared<LinearVelocityIntegrator>(m_strain_parameterisation_stack,
    //                                                m_angular_velocity,
    //                                                m_rod_properties->m_rod_dimensions.m_L)
    // };

    // //  Integrator for the angular accelerations
    // std::shared_ptr<OSNI::ODESolverInterface> m_angular_acceleration {
    //     std::make_shared<AngularAccelerationIntegrator>(m_strain_parameterisation_stack,
    //                                                     m_angular_velocity,
    //                                                     m_rod_properties->m_rod_dimensions.m_L)
    // };

    // //  Integrator for the linear accelerations
    // std::shared_ptr<OSNI::ODESolverInterface> m_linear_acceleration {
    //     std::make_shared<LinearAccelerationIntegrator>(m_strain_parameterisation_stack,
    //                                                    m_angular_velocity,
    //                                                    m_linear_velocity,
    //                                                    m_angular_acceleration,
    //                                                    m_rod_properties->m_rod_dimensions.m_L)
    // };


    // //  Integrator for the internal forces
    // std::shared_ptr<OSNI::ODESolverInterface> m_internal_forces {
    //     std::make_shared<InternalForcesIntegrator>(m_strain_parameterisation_stack,
    //                                                m_rod_properties,
    //                                                m_angular_velocity,
    //                                                m_linear_velocity,
    //                                                m_linear_acceleration,
    //                                                m_quaternion,
    //                                                m_position,
    //                                                m_rod_properties->m_rod_dimensions.m_L)
    // };

    // //  Integrator for the internal couples
    // std::shared_ptr<OSNI::ODESolverInterface> m_internal_couples {
    //     std::make_shared<InternalCouplesIntegrator>(m_strain_parameterisation_stack,
    //                                                 m_rod_properties,
    //                                                 m_angular_velocity,
    //                                                 m_linear_velocity,
    //                                                 m_angular_acceleration,
    //                                                 m_quaternion,
    //                                                 m_position,
    //                                                 m_internal_forces,
    //                                                 m_rod_properties->m_rod_dimensions.m_L)
    // };

    // //  Integrator for the generalised coordinates
    // std::shared_ptr<OSNI::ODESolverInterface> m_generalised_forces {
    //     std::make_shared<GeneralisedForcesIntegrator>(m_strain_parameterisation_stack,
    //                                                   m_internal_couples,
    //                                                   m_internal_forces,
    //                                                   m_rod_properties->m_rod_dimensions.m_L)
    // };




    /*!
     * \brief updateIntegrationDomain this function is used to set the integration upper limit at every integrator
     * \param t_upper_integration_limit is the upper integration limit used to scale the normalised domain [0, 1]
     */
    void updateIntegrationDomain(const double &t_upper_integration_limit=1.0f);
};


}   //  namespace CROSP::numerical_integrators::spectral_method::idm_integrators





#endif // IDM_INTEGRATORS_HPP
