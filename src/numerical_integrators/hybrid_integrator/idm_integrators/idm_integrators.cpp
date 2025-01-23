/**
 * \file idm_integrators.cpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod IDM
 * \date 03-08-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#include "CROSP/numerical_integrators/hybrid_integrator/idm_integrators/idm_integrators.hpp"



namespace CROSP::numerical_integrators::hybrid_integrator::idm_integrators {






// QuaternionIntegrator::QuaternionIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
//                                            const double &t_upper_integration_limit,
//                                            const Eigen::Vector4d &t_initial_condition)
//     : OSNI::ODEA(4,
//                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
//                  t_strain_parameterisation_stack->m_number_of_points,
//                  t_upper_integration_limit),
//       m_strain_parameterisation_stack(t_strain_parameterisation_stack)
// {
//     this->setInitialConditions(t_initial_condition);
// }



// Eigen::MatrixXd QuaternionIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
// {
//     m_A_at_chebychev_point  <<              0             ,   -m_K_stack->at(t_point)(0),   -m_K_stack->at(t_point)(1),   -m_K_stack->at(t_point)(2),
//                                 m_K_stack->at(t_point)(0) ,                0            ,    m_K_stack->at(t_point)(2),   -m_K_stack->at(t_point)(1),
//                                 m_K_stack->at(t_point)(1) ,   -m_K_stack->at(t_point)(2),                0            ,    m_K_stack->at(t_point)(0),
//                                 m_K_stack->at(t_point)(2) ,    m_K_stack->at(t_point)(1),   -m_K_stack->at(t_point)(0),                0            ;

//     return 0.5*m_A_at_chebychev_point;
// }





PositionIntegrator::PositionIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                                       const double &t_upper_integration_limit,
                                       const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEb(3,
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_strain_parameterisation_stack->m_number_of_points,
                 t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_quaternion(t_quaternion_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::VectorXd PositionIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    m_quaternion_at_point = {m_quaternion->getStateAtPoint(t_point)(0),
                            m_quaternion->getStateAtPoint(t_point)(1),
                            m_quaternion->getStateAtPoint(t_point)(2),
                            m_quaternion->getStateAtPoint(t_point)(3)};

    return m_quaternion_at_point.toRotationMatrix() * m_Gamma_stack->at(t_point);
}







AngularVelocityIntegrator::AngularVelocityIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack)
    : m_dot_K_stack( t_strain_parameterisation_stack->m_dot_K_stack.get() )
{}



Eigen::VectorXd AngularVelocityIntegrator::computeParameterVectorAtPoint(const unsigned int t_point)
{
    return m_dot_K_stack->at(t_point);
}





LinearVelocityIntegrator::LinearVelocityIntegrator(const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                                                   std::shared_ptr<const ODECascade> t_angular_velocity_integrator)
    : m_Gamma_stack(t_strain_parameterisation_stack->m_Gamma_stack.get()),
      m_dot_Gamma_stack(t_strain_parameterisation_stack->m_dot_Gamma_stack.get()),
      m_angular_velocity(t_angular_velocity_integrator)
{}





Eigen::VectorXd LinearVelocityIntegrator::computeParameterVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d Gamma = m_Gamma_stack->at(t_point);
    Eigen::Vector3d dot_Gamma = m_dot_Gamma_stack->at(t_point);
    Eigen::Vector3d Omega = m_angular_velocity->getStateAtPoint(t_point);

    return dot_Gamma - Gamma.cross(Omega);
}








IDMIntegrators::IDMIntegrators(const strain_parameterisation_stack::StrainParameterisationStack* t_strain_parameterisation_stack,
                               const rod_properties::RodProperties* t_rod_properties)
    : m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_rod_properties(t_rod_properties)
{

    std::shared_ptr<OSNI::ODECascade> angular_velocity =
        std::make_shared<AngularVelocityIntegrator>(m_strain_parameterisation_stack);


    std::shared_ptr<OSNI::ODECascade> linear_velocity =
        std::make_shared<LinearVelocityIntegrator>(m_strain_parameterisation_stack,
                                                   angular_velocity);



    m_forward_kinematics_integrator.addODE( angular_velocity );
    m_forward_kinematics_integrator.addODE( linear_velocity );

}


// void IDMIntegrators::updateIntegrationDomain(const double &t_upper_integration_limit)
// {
//     m_quaternion->setUpperIntegrationDomain(t_upper_integration_limit);
//     m_position->setUpperIntegrationDomain(t_upper_integration_limit);

//     m_angular_velocity->setUpperIntegrationDomain(t_upper_integration_limit);
//     m_linear_velocity->setUpperIntegrationDomain(t_upper_integration_limit);

//     m_angular_acceleration->setUpperIntegrationDomain(t_upper_integration_limit);
//     m_linear_acceleration->setUpperIntegrationDomain(t_upper_integration_limit);

//     m_internal_forces->setUpperIntegrationDomain(t_upper_integration_limit);
//     m_internal_couples->setUpperIntegrationDomain(t_upper_integration_limit);

//     m_generalised_forces->setUpperIntegrationDomain(t_upper_integration_limit);
// }



}   //  namespace CROSP::numerical_integrators::spectral_method::idm_integrators
