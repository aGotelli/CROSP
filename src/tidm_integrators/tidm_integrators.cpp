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


#include "CROSP/tidm_integrators/tidm_integrators.hpp"

namespace CROSP::tidm_integrators {




DeltaRotation::DeltaRotation(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack)
{}

DeltaRotation::DeltaRotation(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                             const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack)
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







DeltaPosition::DeltaPosition(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_Delta_Lambda_stack(t_strain_parameterisation_Delta->m_Lambda_stack),
      m_delta_rotation(t_delta_rotation)
{}

DeltaPosition::DeltaPosition(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation,
                             const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_Delta_Lambda_stack(t_strain_parameterisation_Delta->m_Lambda_stack),
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
    return -::LieAlgebra::skew( m_Lambda_stack->at(t_point) ) * m_delta_rotation->getStateAtPoint(t_point)
            + m_Delta_Lambda_stack->at(t_point);
}





DeltaAngularVelocity::DeltaAngularVelocity(const unsigned int t_number_of_Chebyshev_points,
                                           std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                           std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                           std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack),
      m_Delta_dot_K_stack(t_strain_parameterisation_Delta->m_dot_K_stack),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity)
{}

DeltaAngularVelocity::DeltaAngularVelocity(const unsigned int t_number_of_Chebyshev_points,
                                           std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                           std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                           std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                           const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack),
      m_Delta_dot_K_stack(t_strain_parameterisation_Delta->m_dot_K_stack),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity)
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


DeltaLinearVelocity::DeltaLinearVelocity(const unsigned int t_number_of_Chebyshev_points,
                                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                         std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                         std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_Delta_Lambda_stack(t_strain_parameterisation_Delta->m_Lambda_stack),
      m_Delta_dot_Lambda_stack(t_strain_parameterisation_Delta->m_dot_Lambda_stack),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity)
{}

DeltaLinearVelocity::DeltaLinearVelocity(const unsigned int t_number_of_Chebyshev_points,
                                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                         std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                         std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                         const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_Delta_Lambda_stack(t_strain_parameterisation_Delta->m_Lambda_stack),
      m_Delta_dot_Lambda_stack(t_strain_parameterisation_Delta->m_dot_Lambda_stack),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaLinearVelocity::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}

Eigen::VectorXd DeltaLinearVelocity::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_Delta_dot_Lambda_stack->at(t_point)
            - ::LieAlgebra::skew( m_Lambda_stack->at(t_point) ) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ) * m_linear_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_Lambda_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point);
}




DeltaAngularAccelerations::DeltaAngularAccelerations(const unsigned int t_number_of_Chebyshev_points,
                                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                                     std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack),
      m_Delta_dot_K_stack(t_strain_parameterisation_Delta->m_dot_K_stack),
      m_Delta_ddot_K_stack(t_strain_parameterisation_Delta->m_ddot_K_stack),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
      m_angular_acceleration_integrator(t_idm_integrators->m_angular_acceleration),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator)
{}


DeltaAngularAccelerations::DeltaAngularAccelerations(const unsigned int t_number_of_Chebyshev_points,
                                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                                     std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                     const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack),
      m_Delta_dot_K_stack(t_strain_parameterisation_Delta->m_dot_K_stack),
      m_Delta_ddot_K_stack(t_strain_parameterisation_Delta->m_ddot_K_stack),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
      m_angular_acceleration_integrator(t_idm_integrators->m_angular_acceleration),
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




DeltaLinearAccelerations::DeltaLinearAccelerations(const unsigned int t_number_of_Chebyshev_points,
                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                                   std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack),
      m_Delta_dot_K_stack(t_strain_parameterisation_Delta->m_dot_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_dot_Lambda_stack(t_strain_parameterisation->m_dot_Lambda_stack),
      m_Delta_Lambda_stack(t_strain_parameterisation_Delta->m_Lambda_stack),
      m_Delta_dot_Lambda_stack(t_strain_parameterisation_Delta->m_dot_Lambda_stack),
      m_Delta_ddot_Lambda_stack(t_strain_parameterisation_Delta->m_ddot_Lambda_stack),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
      m_angular_acceleration_integrator(t_idm_integrators->m_angular_acceleration),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
      m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity),
      m_linear_acceleration_integrator(t_idm_integrators->m_linear_acceleration),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator)
{}

DeltaLinearAccelerations::DeltaLinearAccelerations(const unsigned int t_number_of_Chebyshev_points,
                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                                   std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                   const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_Delta->m_K_stack),
      m_Delta_dot_K_stack(t_strain_parameterisation_Delta->m_dot_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_dot_Lambda_stack(t_strain_parameterisation->m_dot_Lambda_stack),
      m_Delta_Lambda_stack(t_strain_parameterisation_Delta->m_Lambda_stack),
      m_Delta_dot_Lambda_stack(t_strain_parameterisation_Delta->m_dot_Lambda_stack),
      m_Delta_ddot_Lambda_stack(t_strain_parameterisation_Delta->m_ddot_Lambda_stack),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
      m_angular_acceleration_integrator(t_idm_integrators->m_angular_acceleration),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
      m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity),
      m_linear_acceleration_integrator(t_idm_integrators->m_linear_acceleration),
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
    return m_Delta_ddot_Lambda_stack->at(t_point)
            - ::LieAlgebra::skew( m_Lambda_stack->at(t_point) ) * m_Delta_angular_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_Lambda_stack->at(t_point) ) * m_angular_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ) * m_linear_acceleration_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_dot_Lambda_stack->at(t_point) ) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_dot_Lambda_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_Delta_dot_K_stack->at(t_point) ) * m_linear_velocity_integrator->getStateAtPoint(t_point);
}





DeltaInternalForcesIntegrator::DeltaInternalForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                                             std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                             std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_number_of_Chebyshev_points),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_rod_properties(t_rod_properties),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_perturbation->m_K_stack),
      m_quaternion_integrator(t_idm_integrators->m_quaternion),
      m_position_integrator(t_idm_integrators->m_position),
      m_Delta_rotation_integrator(t_Delta_rotation_integrator),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
      m_Delta_linear_acceleration_integrator(t_Delta_linear_acceleration_integrator),
      m_internal_forces_integrator(t_idm_integrators->m_internal_forces)
{}

DeltaInternalForcesIntegrator::DeltaInternalForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                                             std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                             std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator,
                                                             const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_number_of_Chebyshev_points),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_rod_properties(t_rod_properties),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Delta_K_stack(t_strain_parameterisation_perturbation->m_K_stack),
      m_quaternion_integrator(t_idm_integrators->m_quaternion),
      m_position_integrator(t_idm_integrators->m_position),
      m_Delta_rotation_integrator(t_Delta_rotation_integrator),
      m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
      m_Delta_linear_acceleration_integrator(t_Delta_linear_acceleration_integrator),
      m_internal_forces_integrator(t_idm_integrators->m_internal_forces)
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
    if(t_point == 0 or t_point == m_number_of_Chebyshev_points-1)
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







DeltaInternalCouplesIntegrator::DeltaInternalCouplesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                               std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                               std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                                               std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                               std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_number_of_Chebyshev_points),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_rod_properties(t_rod_properties),
      m_idm_integrators(t_idm_integrators),
      m_strain_parameterisation(t_strain_parameterisation),
      m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation),
      m_Delta_rotation_integrator(t_Delta_rotation_integrator),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
      m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator)
{}

DeltaInternalCouplesIntegrator::DeltaInternalCouplesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                               std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                               std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                                               std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                                               std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                                               std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                                               const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_number_of_Chebyshev_points),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_rod_properties(t_rod_properties),
      m_idm_integrators(t_idm_integrators),
      m_strain_parameterisation(t_strain_parameterisation),
      m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation),
      m_Delta_rotation_integrator(t_Delta_rotation_integrator),
      m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
      m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
      m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
      m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::MatrixXd DeltaInternalCouplesIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
}

Eigen::VectorXd DeltaInternalCouplesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{

    return ::LieAlgebra::skew(m_Lambda_stack->at(t_point)).transpose() * m_Delta_internal_forces_integrator->getStateAtPoint(t_point)
            + ::LieAlgebra::skew(m_Delta_K_stack->at(t_point)).transpose() * m_internal_couples_integrator->getStateAtPoint(t_point)
            + ::LieAlgebra::skew(m_Delta_Lambda_stack->at(t_point)).transpose() * m_internal_forces_integrator->getStateAtPoint(t_point)
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




DeltaGeneralisedForcesIntegrator::DeltaGeneralisedForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                                   const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator)
    : OSNI::ODEb(t_polynomial_representation.getCoordinatesDimension(),
                 ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                 t_number_of_Chebyshev_points),
      m_B(t_polynomial_representation.m_B),
      m_strain_parameterisation(t_strain_parameterisation),
      m_Delta_internal_couples_integrator(t_Delta_internal_couples_integrator),
      m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator)
{}


DeltaGeneralisedForcesIntegrator::DeltaGeneralisedForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                                   const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                                                   const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEb(t_polynomial_representation.getCoordinatesDimension(),
                 ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                 t_number_of_Chebyshev_points),
      m_B(t_polynomial_representation.m_B),
      m_strain_parameterisation(t_strain_parameterisation),
      m_Delta_internal_couples_integrator(t_Delta_internal_couples_integrator),
      m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::VectorXd DeltaGeneralisedForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d Delta_C = m_Delta_internal_couples_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d Delta_N = m_Delta_internal_forces_integrator->getStateAtPoint(t_point);

    Eigen::Matrix<double, 6, 1> Delta_Lambda;
    Delta_Lambda << Delta_C,
                    Delta_N;

    const auto Phi = m_strain_parameterisation->m_Phi_stack.at(t_point);
    return -Phi.transpose()*m_B.transpose()*Delta_Lambda;
}





TIDMIntegrators::TIDMIntegrators(const unsigned int t_number_of_Chebyshev_points,
                                 const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                                 std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                 std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                                 std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                 std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_polynomial_representation(t_polynomial_representation),
      m_strain_parameterisation(t_strain_parameterisation),
      m_strain_parameterisation_Delta(t_strain_parameterisation_Delta),
      m_idm_integrators(t_idm_integrators),
      m_rod_properties(t_rod_properties)
{}




}   //  namespace CROSP::tidm_integrators
