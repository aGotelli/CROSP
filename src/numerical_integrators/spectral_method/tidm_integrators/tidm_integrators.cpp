/**
 * \file tidm_integrators.cpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod TIDM
 * \date 03-08-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#include "CROSP/numerical_integrators/spectral_method/tidm_integrators/tidm_integrators.hpp"


namespace CROSP::numerical_integrators::spectral_method::tidm_integrators {




DeltaRotation::DeltaRotation(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                             std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                             const double &t_upper_integration_limit,
                             const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack )
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaRotation::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}

Eigen::VectorXd DeltaRotation::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_Delta_K_stack->at(t_point);
}







DeltaPosition::DeltaPosition(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                             std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation,
                             const double &t_upper_integration_limit,
                             const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack ),
      m_delta_rotation(t_delta_rotation)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaPosition::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}

Eigen::VectorXd DeltaPosition::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return - m_hat_Gamma_stack->at(t_point) * m_delta_rotation->getStateAtPoint(t_point)
            + m_Delta_Gamma_stack->at(t_point);
}




DeltaAngularVelocity::DeltaAngularVelocity(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                           std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                                           std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                           const double &t_upper_integration_limit,
                                           const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack ),
      m_idm_integrators( t_idm_integrators )
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaAngularVelocity::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}

Eigen::VectorXd DeltaAngularVelocity::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return - m_hat_Delta_K_stack->at(t_point)*m_angular_velocity_integrator->getStateAtPoint(t_point)
            + m_Delta_dot_K_stack->at(t_point);
}





DeltaLinearVelocity::DeltaLinearVelocity(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                         std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                                         std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                         std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                         const double &t_upper_integration_limit,
                                         const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack ),
      m_idm_integrators( t_idm_integrators ),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaLinearVelocity::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}

Eigen::VectorXd DeltaLinearVelocity::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_Delta_dot_Gamma_stack->at(t_point)
            - m_hat_Gamma_stack->at(t_point) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - m_hat_Delta_K_stack->at(t_point) * m_linear_velocity_integrator->getStateAtPoint(t_point)
            - m_hat_Delta_Gamma_stack->at(t_point) * m_angular_velocity_integrator->getStateAtPoint(t_point);
}





DeltaAngularAccelerations::DeltaAngularAccelerations(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                     std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                                                     std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                     const double &t_upper_integration_limit,
                                                     const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack ),
      m_idm_integrators( t_idm_integrators ),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaAngularAccelerations::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}

Eigen::VectorXd DeltaAngularAccelerations::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_Delta_ddot_K_stack->at(t_point)
            - m_hat_Delta_K_stack->at(t_point) * m_angular_acceleration_integrator->getStateAtPoint(t_point)
            - m_hat_dot_K_stack->at(t_point) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - m_hat_Delta_dot_K_stack->at(t_point) * m_angular_velocity_integrator->getStateAtPoint(t_point);
}





DeltaLinearAccelerations::DeltaLinearAccelerations(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                   std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                                                   std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                   const double &t_upper_integration_limit,
                                                   const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack ),
      m_idm_integrators( t_idm_integrators ),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaLinearAccelerations::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}

Eigen::VectorXd DeltaLinearAccelerations::computerParametersVectorAtPoint(const unsigned int t_point)
{
//    return m_Delta_ddot_Gamma_stack->at(t_point)
//            - m_hat_Gamma_stack->at(t_point) * m_Delta_angular_acceleration_integrator->getStateAtPoint(t_point)
//            - m_hat_Delta_Gamma_stack->at(t_point) * m_angular_acceleration_integrator->getStateAtPoint(t_point)
//            - m_hat_Delta_K_stack->at(t_point) * m_linear_acceleration_integrator->getStateAtPoint(t_point)
//            - m_hat_dot_Gamma_stack->at(t_point) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
//            - m_hat_dot_K_stack->at(t_point) * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
//            - m_hat_Delta_dot_Gamma_stack->at(t_point) * m_angular_velocity_integrator->getStateAtPoint(t_point)
//            - m_hat_Delta_dot_K_stack->at(t_point) * m_linear_velocity_integrator->getStateAtPoint(t_point);

    Eigen::Vector3d Gamma_Delta_dot_Omega = m_hat_Gamma_stack->at(t_point) * m_Delta_angular_acceleration_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d dot_Gamma_Delta_Omega = m_hat_dot_Gamma_stack->at(t_point) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d dot_K_Delta_V = m_hat_dot_K_stack->at(t_point) * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d V_Delta_dot_K = - m_hat_Delta_dot_K_stack->at(t_point) * m_linear_velocity_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d Omega_Delta_dot_Gamma = - m_hat_Delta_dot_Gamma_stack->at(t_point) * m_angular_velocity_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d dot_V_Delta_K = - m_hat_Delta_K_stack->at(t_point) * m_linear_acceleration_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d dot_Omega_Delta_Gamma = - m_hat_Delta_dot_Gamma_stack->at(t_point) * m_angular_acceleration_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d ddot_Gamma = m_Delta_ddot_Gamma_stack->at(t_point);

//    if(t_point == 15)
//        std::cout << "Point : " << t_point << "\n"
//                     "  Gamma_Delta_dot_Omega : \n" << Gamma_Delta_dot_Omega << "\n"
//                     "  dot_Gamma_Delta_Omega : \n" << dot_Gamma_Delta_Omega << "\n"
//                     "  dot_K_Delta_V : \n" << dot_K_Delta_V << "\n"
//                     "  V_Delta_dot_K : \n" << V_Delta_dot_K << "\n"
//                     "  Omega_Delta_dot_Gamma : \n" << Omega_Delta_dot_Gamma << "\n"
//                     "  dot_V_Delta_K : \n" << dot_V_Delta_K << "\n"
//                     "  dot_Omega_Delta_Gamma : \n" << dot_Omega_Delta_Gamma << "\n"
//                     "  ddot_Gamma : \n" << ddot_Gamma << "\n"
//                     "\n\n";


    Eigen::Vector3d b = - Gamma_Delta_dot_Omega
            - dot_Gamma_Delta_Omega
            - dot_K_Delta_V
            + V_Delta_dot_K
            + Omega_Delta_dot_Gamma
            + dot_V_Delta_K
            + dot_Omega_Delta_Gamma
            + ddot_Gamma;

    return b;
}





DeltaInternalForcesIntegrator::DeltaInternalForcesIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                             std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                                                             std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                             std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator,
                                                             const double &t_upper_integration_limit,
                                                             const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack ),
      m_idm_integrators( t_idm_integrators ),
      m_Delta_rotation_integrator(t_Delta_rotation_integrator),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
      m_Delta_linear_acceleration_integrator(t_Delta_linear_acceleration_integrator),
      m_rod_properties(t_rod_properties)
{
 this->setInitialConditions(t_initial_condition);
}




Eigen::MatrixXd DeltaInternalForcesIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return m_hat_K_stack->at(t_point).transpose();
}

Eigen::VectorXd DeltaInternalForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    const auto M_linear = m_rod_properties->getMLinear();

    return M_linear*m_Delta_linear_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * M_linear * m_linear_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * M_linear * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
            - computeLocalExternalForces(t_point)
            +  m_hat_Delta_K_stack->at(t_point).transpose()*m_internal_forces_integrator->getStateAtPoint(t_point);

}


Eigen::Vector3d DeltaInternalForcesIntegrator::computeLocalExternalForces(unsigned int t_point) const
{
    if(t_point == 0 or t_point == m_strain_parameterisation_stack->m_number_of_points-1)
        return Eigen::Vector3d::Zero();

    Eigen::Vector4d q(m_quaternion_integrator->getStateAtPoint(t_point));

    Eigen::Matrix3d R = Eigen::Quaterniond(q(0),
                                           q(1),
                                           q(2),
                                           q(3)).toRotationMatrix();

    Eigen::Vector3d Delta_rotation = m_Delta_rotation_integrator->getStateAtPoint(t_point);

    Eigen::Vector3d Delta_N_bar = ::LieAlgebra::skew( Delta_rotation ).transpose()
                                    *R.transpose()
                                    *m_rod_properties->distributedGravitationalForce();


    return Delta_N_bar;
}






DeltaInternalCouplesIntegrator::DeltaInternalCouplesIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                               std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                                                               std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                               std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                                               const double &t_upper_integration_limit,
                                                               const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack ),
      m_idm_integrators(t_idm_integrators),
      m_Delta_rotation_integrator(t_Delta_rotation_integrator),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
      m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator),
      m_rod_properties(t_rod_properties)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaInternalCouplesIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return m_hat_K_stack->at(t_point).transpose();
}

Eigen::VectorXd DeltaInternalCouplesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{

    const auto M_linear = m_rod_properties->getMLinear();
    const auto M_angular = m_rod_properties->getMAngular();

    return m_hat_Gamma_stack->at(t_point).transpose() * m_Delta_internal_forces_integrator->getStateAtPoint(t_point)
            + m_hat_Delta_K_stack->at(t_point).transpose() * m_internal_couples_integrator->getStateAtPoint(t_point)
            + m_hat_Delta_Gamma_stack->at(t_point).transpose() * m_internal_forces_integrator->getStateAtPoint(t_point)
            + M_angular * m_Delta_angular_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_angular_velocity_integrator->getStateAtPoint(t_point)).transpose() * M_angular * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_linear_velocity_integrator->getStateAtPoint(t_point)).transpose() * M_linear * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)).transpose() * M_angular * m_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)).transpose() * M_linear * m_linear_velocity_integrator->getStateAtPoint(t_point)
            - computeLocalExternalCouples(t_point);

}


Eigen::Vector3d DeltaInternalCouplesIntegrator::computeLocalExternalCouples(unsigned int) const
{
    return Eigen::Vector3d::Zero();
}







DeltaGeneralisedForcesIntegrator::DeltaGeneralisedForcesIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                                                   const double &t_upper_integration_limit)
    : OSNI::ODEb(t_strain_parameterisation_stack->m_B_Phi_stack.cols(),
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_strain_parameterisation_stack->m_number_of_points,
                 t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_coordinates_dimension(t_strain_parameterisation_stack->m_B_Phi_stack.cols()),
      m_Delta_internal_couples_integrator(t_Delta_internal_couples_integrator),
      m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator)
{
    this->setInitialConditions( Eigen::VectorXd::Zero(m_coordinates_dimension) );
}

Eigen::VectorXd DeltaGeneralisedForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d Delta_C = m_Delta_internal_couples_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d Delta_N = m_Delta_internal_forces_integrator->getStateAtPoint(t_point);

    Eigen::Matrix<double, 6, 1> Delta_Lambda;
    Delta_Lambda << Delta_C,
                    Delta_N;


    const auto b_at_point = m_strain_parameterisation_stack->m_B_Phi_stack.block(6*t_point, 0, 6, m_coordinates_dimension).transpose()*Delta_Lambda;
    return b_at_point;
}





TIDMIntegrators::TIDMIntegrators(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                 std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                                 std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                 std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
    : m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_Delta_strain_parameterisation_stack( t_Delta_strain_parameterisation_stack ),
      m_idm_integrators(t_idm_integrators),
      m_rod_properties(t_rod_properties)
{}


void TIDMIntegrators::updateIntegrationDomain(const double &t_upper_integration_limit)
    {
        m_Delta_rotation->setUpperIntegrationDomain(t_upper_integration_limit);
        m_Delta_position->setUpperIntegrationDomain(t_upper_integration_limit);

        m_Delta_angular_velocity->setUpperIntegrationDomain(t_upper_integration_limit);
        m_Delta_linear_velocity->setUpperIntegrationDomain(t_upper_integration_limit);

        m_Delta_angular_acceleration->setUpperIntegrationDomain(t_upper_integration_limit);
        m_Delta_linear_acceleration->setUpperIntegrationDomain(t_upper_integration_limit);

        m_Delta_internal_forces->setUpperIntegrationDomain(t_upper_integration_limit);
        m_Delta_internal_couples->setUpperIntegrationDomain(t_upper_integration_limit);

        m_Delta_generalised_forces->setUpperIntegrationDomain(t_upper_integration_limit);
    }



}   //  namespace CROSP::numerical_integrators::spectral_method::tidm_integrators
