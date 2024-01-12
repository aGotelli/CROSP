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
#include "ATORS/ATORS.hpp"
#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"



/// \brief CROSP::idm_integrators is the namespace containing the integrators for the IDM
namespace CROSP::numerical_integrators::spectral_method::idm_integrators {







/*!
 * \brief The QuaternionIntegrator struct implements the spectral numerical integration for the quaternion
 *
 * This class implemnts the integration for the quaternions forllowing the formula \f$ Q' = \frac{1}{2}  A_{(K)} Q
 */
struct QuaternionIntegrator : public OSNI::ODEA {

    QuaternionIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                         const double &t_upper_integration_limit=1.0f,
                         const Eigen::Vector4d &t_initial_condition=Eigen::Vector4d(1, 0, 0, 0));

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;



    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_K_stack {
        m_strain_parameterisation_stack->m_K_stack.get()
    };

    Eigen::Matrix4d m_A_at_chebychev_point;

};



struct PositionIntegrator : public OSNI::ODEb {

    PositionIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                       const double &t_upper_integration_limit=1.0f,
                       const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_Gamma_stack {
        m_strain_parameterisation_stack->m_Gamma_stack.get()
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    Eigen::Quaterniond m_quaternion_at_point;

};




struct AngularVelocityIntegrator : public OSNI::ODEAb {

    AngularVelocityIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                              const double &t_upper_integration_limit=1.0f,
                              const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_K_stack {
        m_strain_parameterisation_stack->m_K_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_dot_K_stack {
        m_strain_parameterisation_stack->m_dot_K_stack.get()
    };


};


struct LinearVelocityIntegrator : public OSNI::ODEAb {

    LinearVelocityIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             const double &t_upper_integration_limit=1.0f,
                             const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point);


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point);


    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_K_stack {
        m_strain_parameterisation_stack->m_K_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_Gamma_stack {
        m_strain_parameterisation_stack->m_Gamma_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_dot_Gamma_stack {
        m_strain_parameterisation_stack->m_dot_Gamma_stack.get()
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;



};




struct AngularAccelerationIntegrator : public OSNI::ODEAb {

    AngularAccelerationIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                  const double &t_upper_integration_limit=1.0f,
                                  const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_K_stack {
        m_strain_parameterisation_stack->m_K_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_dot_K_stack {
        m_strain_parameterisation_stack->m_dot_K_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_ddot_K_stack {
        m_strain_parameterisation_stack->m_ddot_K_stack.get()
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;


};







struct LinearAccelerationIntegrator : public OSNI::ODEAb {

    LinearAccelerationIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                 const double &t_upper_integration_limit=1.0f,
                                 const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_K_stack {
        m_strain_parameterisation_stack->m_K_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_dot_K_stack {
        m_strain_parameterisation_stack->m_dot_K_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_Gamma_stack {
        m_strain_parameterisation_stack->m_Gamma_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_dot_Gamma_stack{
        m_strain_parameterisation_stack->m_dot_Gamma_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_ddot_Gamma_stack{
        m_strain_parameterisation_stack->m_ddot_Gamma_stack.get()
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;



};




struct InternalForcesIntegrator : public OSNI::ODEAb {

    InternalForcesIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                             const rod_properties::RodProperties* t_rod_properties,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                             const double &t_upper_integration_limit=1.0f,
                             const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computeDistributedForce(const unsigned int t_point) const;




    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_K_stack {
        m_strain_parameterisation_stack->m_K_stack.get()
    };

    const rod_properties::RodProperties* m_rod_properties;

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position;

};


struct InternalCouplesIntegrator : public OSNI::ODEAb {

    InternalCouplesIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                              const rod_properties::RodProperties* t_rod_properties,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                              const double &t_upper_integration_limit=1.0f,
                              const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());




    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computeDistributedCouple(const unsigned int) const;


    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_K_stack {
        m_strain_parameterisation_stack->m_K_stack.get()
    };

    const ::CROSP::strain_parameterisation_stack::EigenRef3Stack* m_Gamma_stack {
        m_strain_parameterisation_stack->m_Gamma_stack.get()
    };

    const rod_properties::RodProperties* m_rod_properties;

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position;

    std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces;

};


struct GeneralisedForcesIntegrator : public OSNI::ODEb {

    GeneralisedForcesIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                                std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
                                std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                                const double &t_upper_integration_limit=1.0f);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* m_strain_parameterisation_stack;

    unsigned int m_coordinates_dimension;

    std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator;
    std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;

};







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


    //  Integrator for the quaternions
    std::shared_ptr<OSNI::ODESolverInterface> m_quaternion {
        std::make_shared<QuaternionIntegrator>(m_strain_parameterisation_stack,
                                               m_rod_properties->m_rod_dimensions.m_L)
    };

    //  Integrator for the positions
    std::shared_ptr<OSNI::ODESolverInterface> m_position {
        std::make_shared<PositionIntegrator>(m_strain_parameterisation_stack,
                                             m_quaternion,
                                             m_rod_properties->m_rod_dimensions.m_L)
    };

    //  Integrator for the angular velocities
    std::shared_ptr<OSNI::ODESolverInterface> m_angular_velocity {
        std::make_shared<AngularVelocityIntegrator>(m_strain_parameterisation_stack,
                                                    m_rod_properties->m_rod_dimensions.m_L)
    };

    //  Integrator for the linear velocities
    std::shared_ptr<OSNI::ODESolverInterface> m_linear_velocity {
        std::make_shared<LinearVelocityIntegrator>(m_strain_parameterisation_stack,
                                                   m_angular_velocity,
                                                   m_rod_properties->m_rod_dimensions.m_L)
    };

    //  Integrator for the angular accelerations
    std::shared_ptr<OSNI::ODESolverInterface> m_angular_acceleration {
        std::make_shared<AngularAccelerationIntegrator>(m_strain_parameterisation_stack,
                                                        m_angular_velocity,
                                                        m_rod_properties->m_rod_dimensions.m_L)
    };

    //  Integrator for the linear accelerations
    std::shared_ptr<OSNI::ODESolverInterface> m_linear_acceleration {
        std::make_shared<LinearAccelerationIntegrator>(m_strain_parameterisation_stack,
                                                       m_angular_velocity,
                                                       m_linear_velocity,
                                                       m_angular_acceleration,
                                                       m_rod_properties->m_rod_dimensions.m_L)
    };


    //  Integrator for the internal forces
    std::shared_ptr<OSNI::ODESolverInterface> m_internal_forces {
        std::make_shared<InternalForcesIntegrator>(m_strain_parameterisation_stack,
                                                   m_rod_properties,
                                                   m_angular_velocity,
                                                   m_linear_velocity,
                                                   m_linear_acceleration,
                                                   m_quaternion,
                                                   m_position,
                                                   m_rod_properties->m_rod_dimensions.m_L)
    };

    //  Integrator for the internal couples
    std::shared_ptr<OSNI::ODESolverInterface> m_internal_couples {
        std::make_shared<InternalCouplesIntegrator>(m_strain_parameterisation_stack,
                                                    m_rod_properties,
                                                    m_angular_velocity,
                                                    m_linear_velocity,
                                                    m_angular_acceleration,
                                                    m_quaternion,
                                                    m_position,
                                                    m_internal_forces,
                                                    m_rod_properties->m_rod_dimensions.m_L)
    };

    //  Integrator for the generalised coordinates
    std::shared_ptr<OSNI::ODESolverInterface> m_generalised_forces {
        std::make_shared<GeneralisedForcesIntegrator>(m_strain_parameterisation_stack,
                                                      m_internal_couples,
                                                      m_internal_forces,
                                                      m_rod_properties->m_rod_dimensions.m_L)
    };




    /*!
     * \brief updateIntegrationDomain this function is used to set the integration upper limit at every integrator
     * \param t_upper_integration_limit is the upper integration limit used to scale the normalised domain [0, 1]
     */
    void updateIntegrationDomain(const double &t_upper_integration_limit=1.0f);
};


}   //  namespace CROSP::numerical_integrators::spectral_method::idm_integrators





#endif // IDM_INTEGRATORS_HPP
