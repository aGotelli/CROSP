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


#include "CROSP/numerical_integrators/spectral/tidm_integrators.hpp"

namespace CROSP::spectral_integrators::tidm_integrators {





DeltaRotation::DeltaRotation(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                             std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                             const double &t_upper_integration_limit,
                             const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_parameterisation_stack->m_number_of_Chebyshev_points,
                  t_upper_integration_limit),
      m_parameterisation_stack( t_parameterisation_stack ),
      m_parameterisation_stack_Delta( t_parameterisation_stack_Delta )
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaRotation::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}

Eigen::VectorXd DeltaRotation::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_Delta_K_stack->at(t_point);
}







DeltaPosition::DeltaPosition(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                             std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation,
                             const double &t_upper_integration_limit,
                             const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_parameterisation_stack->m_number_of_Chebyshev_points,
                  t_upper_integration_limit),
      m_parameterisation_stack( t_parameterisation_stack ),
      m_parameterisation_stack_Delta( t_parameterisation_stack_Delta ),
      m_delta_rotation(t_delta_rotation)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaPosition::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}

Eigen::VectorXd DeltaPosition::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_Gamma_stack->at(t_point) ) * m_delta_rotation->getStateAtPoint(t_point)
            + m_Delta_Gamma_stack->at(t_point);
}




DeltaAngularVelocity::DeltaAngularVelocity(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                           std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                           std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                           const double &t_upper_integration_limit,
                                           const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_parameterisation_stack->m_number_of_Chebyshev_points,
                  t_upper_integration_limit),
      m_parameterisation_stack( t_parameterisation_stack ),
      m_parameterisation_stack_Delta( t_parameterisation_stack_Delta ),
      m_idm_integrators( t_idm_integrators )
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaAngularVelocity::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}

Eigen::VectorXd DeltaAngularVelocity::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_Delta_K_stack->at(t_point) )*m_angular_velocity_integrator->getStateAtPoint(t_point)
     +m_Delta_dot_K_stack->at(t_point);
}





DeltaLinearVelocity::DeltaLinearVelocity(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                         std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                         std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                         std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                         const double &t_upper_integration_limit,
                                         const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_parameterisation_stack->m_number_of_Chebyshev_points,
                  t_upper_integration_limit),
      m_parameterisation_stack( t_parameterisation_stack ),
      m_parameterisation_stack_Delta( t_parameterisation_stack_Delta ),
      m_idm_integrators( t_idm_integrators ),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaLinearVelocity::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}

Eigen::VectorXd DeltaLinearVelocity::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_Delta_dot_Gamma_stack->at(t_point)
            - ::LieAlgebra::skew( m_Gamma_stack->at(t_point) ) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ) * m_linear_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_Gamma_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point);
}





DeltaAngularAccelerations::DeltaAngularAccelerations(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                                     std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                                     std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                     const double &t_upper_integration_limit,
                                                     const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_parameterisation_stack->m_number_of_Chebyshev_points,
                  t_upper_integration_limit),
      m_parameterisation_stack( t_parameterisation_stack ),
      m_parameterisation_stack_Delta( t_parameterisation_stack_Delta ),
      m_idm_integrators( t_idm_integrators ),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaAngularAccelerations::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}

Eigen::VectorXd DeltaAngularAccelerations::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_Delta_ddot_K_stack->at(t_point)
            - ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ) * m_angular_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_dot_K_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point);
}





DeltaLinearAccelerations::DeltaLinearAccelerations(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                                   std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                                   std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                   const double &t_upper_integration_limit,
                                                   const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_parameterisation_stack->m_number_of_Chebyshev_points,
                  t_upper_integration_limit),
      m_parameterisation_stack( t_parameterisation_stack ),
      m_parameterisation_stack_Delta( t_parameterisation_stack_Delta ),
      m_idm_integrators( t_idm_integrators ),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaLinearAccelerations::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}

Eigen::VectorXd DeltaLinearAccelerations::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_Delta_ddot_Gamma_stack->at(t_point)
            - ::LieAlgebra::skew( m_Gamma_stack->at(t_point) ) * m_Delta_angular_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_Gamma_stack->at(t_point) ) * m_angular_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ) * m_linear_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_dot_Gamma_stack->at(t_point) ) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_dot_Gamma_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_dot_K_stack->at(t_point) ) * m_linear_velocity_integrator->getStateAtPoint(t_point);
}





DeltaInternalForcesIntegrator::DeltaInternalForcesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                                             std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                                             std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                             const rod_properties::RodProperties t_rod_properties,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator,
                                                             const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_parameterisation_stack->m_number_of_Chebyshev_points,
                  t_rod_properties.m_rod_dimensions.m_L),
      m_parameterisation_stack( t_parameterisation_stack ),
      m_parameterisation_stack_Delta( t_parameterisation_stack_Delta ),
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
    return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
}

Eigen::VectorXd DeltaInternalForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{

    return m_M_linear*m_Delta_linear_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_M_linear * m_linear_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_M_linear * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
            - computeLocalExternalForces(t_point)
            + ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ).transpose()*m_internal_forces_integrator->getStateAtPoint(t_point);

}


Eigen::Vector3d DeltaInternalForcesIntegrator::computeLocalExternalForces(unsigned int t_point) const
{
    if(t_point == 0 or t_point == m_parameterisation_stack->m_number_of_Chebyshev_points-1)
        return Eigen::Vector3d::Zero();

    Eigen::Vector4d q(m_quaternion_integrator->getStateAtPoint(t_point));

    Eigen::Matrix3d R = Eigen::Quaterniond(q(0),
                                           q(1),
                                           q(2),
                                           q(3)).toRotationMatrix();

    Eigen::Vector3d Delta_rotation = m_Delta_rotation_integrator->getStateAtPoint(t_point);

    Eigen::Vector3d Delta_N_bar = ::LieAlgebra::skew( Delta_rotation ).transpose()
                                    *R.transpose()
                                    *m_rod_properties.distributedGravitationalForce();


    return Delta_N_bar;
}






DeltaInternalCouplesIntegrator::DeltaInternalCouplesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                                               std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                                               std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                               const rod_properties::RodProperties t_rod_properties,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                                               const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_parameterisation_stack->m_number_of_Chebyshev_points,
                  t_rod_properties.m_rod_dimensions.m_L),
      m_parameterisation_stack( t_parameterisation_stack ),
      m_parameterisation_stack_Delta( t_parameterisation_stack_Delta ),
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
    return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
}

Eigen::VectorXd DeltaInternalCouplesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{

    return ::LieAlgebra::skew(m_Gamma_stack->at(t_point)).transpose() * m_Delta_internal_forces_integrator->getStateAtPoint(t_point)
            + ::LieAlgebra::skew(m_Delta_K_stack->at(t_point)).transpose() * m_internal_couples_integrator->getStateAtPoint(t_point)
            + ::LieAlgebra::skew(m_Delta_Gamma_stack->at(t_point)).transpose() * m_internal_forces_integrator->getStateAtPoint(t_point)
            + m_M_angular * m_Delta_angular_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_angular_velocity_integrator->getStateAtPoint(t_point)).transpose() * m_M_angular * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_linear_velocity_integrator->getStateAtPoint(t_point)).transpose() * m_M_linear * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)).transpose() * m_M_angular * m_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)).transpose() * m_M_linear * m_linear_velocity_integrator->getStateAtPoint(t_point)
            - computeLocalExternalCouples(t_point);

}


Eigen::Vector3d DeltaInternalCouplesIntegrator::computeLocalExternalCouples(unsigned int) const
{
    return Eigen::Vector3d::Zero();
}







DeltaGeneralisedForcesIntegrator::DeltaGeneralisedForcesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                                                   const double &t_upper_integration_limit)
    : OSNI::ODEb(t_parameterisation_stack->m_polynomial_representation.getCoordinatesDimension(),
                 ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                 t_parameterisation_stack->m_number_of_Chebyshev_points,
                 t_upper_integration_limit),
      m_parameterisation_stack(t_parameterisation_stack),
      m_Delta_internal_couples_integrator(t_Delta_internal_couples_integrator),
      m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator)
{
    this->setInitialConditions( Eigen::VectorXd::Zero(t_parameterisation_stack->m_polynomial_representation.getCoordinatesDimension()) );
}

Eigen::VectorXd DeltaGeneralisedForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d Delta_C = m_Delta_internal_couples_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d Delta_N = m_Delta_internal_forces_integrator->getStateAtPoint(t_point);

    Eigen::Matrix<double, 6, 1> Delta_Lambda;
    Delta_Lambda << Delta_C,
                    Delta_N;


    const auto b_at_point = -m_parameterisation_stack->m_map_to_strain_stack.at(t_point).transpose()*Delta_Lambda;
    return b_at_point;
}




TIDMIntegrators::TIDMIntegrators(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                 std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                 const std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                 const rod_properties::RodProperties t_rod_properties)
    : m_rod_properties(t_rod_properties),
      m_parameterisation_stack(t_parameterisation_stack),
      m_parameterisation_stack_Delta(t_parameterisation_stack_Delta),
      m_idm_integrators(t_idm_integrators)
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



}   //  namespace CROSP::spectral_integrators::tidm_integrators
