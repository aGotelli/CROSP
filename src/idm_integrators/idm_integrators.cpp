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


#include "CROSP/idm_integrators/idm_integrators.hpp"



namespace CROSP::idm_integrators {



QuaternionIntegrator::QuaternionIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                           std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation)
    : OSNI::ODEA(4,
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack)
{}


QuaternionIntegrator::QuaternionIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                           std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                           const Eigen::Vector4d t_initial_condition)
    : OSNI::ODEA(4,
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack)
{
    this->setInitialConditions(t_initial_condition);
}



Eigen::MatrixXd QuaternionIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    m_A_at_chebychev_point  <<              0             ,   -m_K_stack->at(t_point)(0),   -m_K_stack->at(t_point)(1),   -m_K_stack->at(t_point)(2),
                                m_K_stack->at(t_point)(0) ,                0            ,    m_K_stack->at(t_point)(2),   -m_K_stack->at(t_point)(1),
                                m_K_stack->at(t_point)(1) ,   -m_K_stack->at(t_point)(2),                0            ,    m_K_stack->at(t_point)(0),
                                m_K_stack->at(t_point)(2) ,    m_K_stack->at(t_point)(1),   -m_K_stack->at(t_point)(0),                0            ;

    return 0.5*m_A_at_chebychev_point;
}





PositionIntegrator::PositionIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                       std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator)
    : OSNI::ODEb(3,
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_number_of_Chebyshev_points),
      m_quaternion(t_quaternion_integrator),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack)
{}

PositionIntegrator::PositionIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                       std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                                       const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEb(3,
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_number_of_Chebyshev_points),
      m_quaternion(t_quaternion_integrator),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::VectorXd PositionIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    m_quaternion_at_point = {m_quaternion->getStateAtPoint(t_point)(0),
                            m_quaternion->getStateAtPoint(t_point)(1),
                            m_quaternion->getStateAtPoint(t_point)(2),
                            m_quaternion->getStateAtPoint(t_point)(3)};

    return m_quaternion_at_point.toRotationMatrix() * m_Lambda_stack->at(t_point);
}






AngularVelocityIntegrator::AngularVelocityIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack)
{}

AngularVelocityIntegrator::AngularVelocityIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                     const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack)
{
    this->setInitialConditions(t_initial_condition);
}


Eigen::MatrixXd AngularVelocityIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}


Eigen::VectorXd AngularVelocityIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_dot_K_stack->at(t_point);
}




LinearVelocityIntegrator::LinearVelocityIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_dot_Lambda_stack(t_strain_parameterisation->m_dot_Lambda_stack),
      m_angular_velocity(t_angular_velocity_integrator)
{}

LinearVelocityIntegrator::LinearVelocityIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                   const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_dot_Lambda_stack(t_strain_parameterisation->m_dot_Lambda_stack),
      m_angular_velocity(t_angular_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}


Eigen::MatrixXd LinearVelocityIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}


Eigen::VectorXd LinearVelocityIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d dot_Lambda = m_dot_Lambda_stack->at(t_point);
    Eigen::Matrix3d skew_Lambda = ::LieAlgebra::skew( m_Lambda_stack->at(t_point) );
    Eigen::Vector3d Omega = m_angular_velocity->getStateAtPoint(t_point);

    return dot_Lambda - skew_Lambda*Omega;
}





AngularAccelerationIntegrator::AngularAccelerationIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_ddot_K_stack(t_strain_parameterisation->m_ddot_K_stack),
      m_angular_velocity(t_angular_velocity_integrator)
{}

AngularAccelerationIntegrator::AngularAccelerationIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                             const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_ddot_K_stack(t_strain_parameterisation->m_ddot_K_stack),
      m_angular_velocity(t_angular_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}



Eigen::MatrixXd AngularAccelerationIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}


Eigen::VectorXd AngularAccelerationIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_ddot_K_stack->at(t_point)
            - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_angular_velocity->getStateAtPoint(t_point);
}









LinearAccelerationIntegrator::LinearAccelerationIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                           std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_dot_Lambda_stack(t_strain_parameterisation->m_dot_Lambda_stack),
      m_ddot_Lambda_stack(t_strain_parameterisation->m_ddot_Lambda_stack),
      m_angular_velocity(t_angular_velocity_integrator),
      m_linear_velocity(t_linear_velocity_integrator),
      m_angular_acceleration(t_angular_acceleration_integrator)
{}

LinearAccelerationIntegrator::LinearAccelerationIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                           std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                                           const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_number_of_Chebyshev_points),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_dot_Lambda_stack(t_strain_parameterisation->m_dot_Lambda_stack),
      m_ddot_Lambda_stack(t_strain_parameterisation->m_ddot_Lambda_stack),
      m_angular_velocity(t_angular_velocity_integrator),
      m_linear_velocity(t_linear_velocity_integrator),
      m_angular_acceleration(t_angular_acceleration_integrator)
{
    this->setInitialConditions(t_initial_condition);
}



Eigen::MatrixXd LinearAccelerationIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return -::LieAlgebra::skew( m_K_stack->at(t_point) );
}


Eigen::VectorXd LinearAccelerationIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_ddot_Lambda_stack->at(t_point)
            - ::LieAlgebra::skew( m_Lambda_stack->at(t_point) ) * m_angular_acceleration->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_dot_Lambda_stack->at(t_point) ) * m_angular_velocity->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_linear_velocity->getStateAtPoint(t_point);
}




InternalForcesIntegrator::InternalForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                   std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_number_of_Chebyshev_points),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_rod_properties(t_rod_properties),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_angular_velocity(t_angular_velocity_integrator),
      m_linear_velocity(t_linear_velocity_integrator),
      m_linear_acceleration(t_linear_acceleration_integrator),
      m_quaternion(t_quaternion_integrator),
      m_position(t_position_integrator)
{}

InternalForcesIntegrator::InternalForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                   std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                                                   const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_number_of_Chebyshev_points),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_rod_properties(t_rod_properties),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_angular_velocity(t_angular_velocity_integrator),
      m_linear_velocity(t_linear_velocity_integrator),
      m_linear_acceleration(t_linear_acceleration_integrator),
      m_quaternion(t_quaternion_integrator),
      m_position(t_position_integrator)
{
    this->setInitialConditions(t_initial_condition);
}


Eigen::MatrixXd InternalForcesIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
}


Eigen::VectorXd InternalForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d inertial_acceleration = m_M_linear*m_linear_acceleration->getStateAtPoint(t_point);
    Eigen::Vector3d inertial_velocities = ::LieAlgebra::skew( m_angular_velocity->getStateAtPoint(t_point) ).transpose() * m_M_linear * m_linear_velocity->getStateAtPoint(t_point);
    Eigen::Vector3d distributed_forces = computeDistributedForce(t_point);

    Eigen::Vector3d b = inertial_acceleration
                        - inertial_velocities
                        - distributed_forces
                        ;


    return b;
}


Eigen::VectorXd InternalForcesIntegrator::computeDistributedForce(const unsigned int t_point) const
{
    if(t_point == 0 or t_point == m_number_of_Chebyshev_points-1)
        return Eigen::Vector3d::Zero();

    Eigen::Vector4d q(m_quaternion->getStateAtPoint(t_point));

    Eigen::Matrix3d R = Eigen::Quaterniond(q(0),
                                           q(1),
                                           q(2),
                                           q(3)).toRotationMatrix();


    const Eigen::Vector3d distributed_weight_force = -m_rod_properties->distributedDensity()*Eigen::Vector3d(0, 0, 9.81);
    Eigen::Vector3d N_bar = R.transpose()*distributed_weight_force;

    return N_bar;
}



InternalCouplesIntegrator::InternalCouplesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                     std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_number_of_Chebyshev_points),
      m_rod_properties( t_rod_properties ),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_angular_velocity(t_angular_velocity_integrator),
      m_linear_velocity(t_linear_velocity_integrator),
      m_angular_acceleration(t_angular_acceleration_integrator),
      m_quaternion(t_quaternion_integrator),
      m_position(t_position_integrator),
      m_internal_forces(t_internal_forces_integrator)
{}

InternalCouplesIntegrator::InternalCouplesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                     std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                                                     const Eigen::Vector3d t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_number_of_Chebyshev_points),
      m_rod_properties( t_rod_properties ),
      m_K_stack(t_strain_parameterisation->m_K_stack),
      m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
      m_angular_velocity(t_angular_velocity_integrator),
      m_linear_velocity(t_linear_velocity_integrator),
      m_angular_acceleration(t_angular_acceleration_integrator),
      m_quaternion(t_quaternion_integrator),
      m_position(t_position_integrator),
      m_internal_forces(t_internal_forces_integrator)
{
    this->setInitialConditions(t_initial_condition);
}




Eigen::MatrixXd InternalCouplesIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
}


Eigen::VectorXd InternalCouplesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{

    return ::LieAlgebra::skew( m_Lambda_stack->at(t_point) ).transpose()*m_internal_forces->getStateAtPoint(t_point)
            + m_M_angular*m_angular_acceleration->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_angular_velocity->getStateAtPoint(t_point)).transpose() * m_M_angular * m_angular_velocity->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_linear_velocity->getStateAtPoint(t_point) ).transpose()* m_M_linear * m_linear_velocity->getStateAtPoint(t_point)
            + computeDistributedCouple(t_point);
}


Eigen::VectorXd InternalCouplesIntegrator::computeDistributedCouple(const unsigned int) const
{
    return Eigen::Vector3d::Zero();
}



GeneralisedForcesIntegrator::GeneralisedForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                         const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                                                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                         std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
                                                         std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator)
    : OSNI::ODEb(t_polynomial_representation.getCoordinatesDimension(),
                 ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                 t_number_of_Chebyshev_points),
      m_B(t_polynomial_representation.m_B),
      m_strain_parameterisation(t_strain_parameterisation),
      m_internal_couples_integrator(t_internal_couples_integrator),
      m_internal_forces_integrator(t_internal_forces_integrator)
{}

GeneralisedForcesIntegrator::GeneralisedForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                                         const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                                                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                         std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
                                                         std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                                                         const Eigen::VectorXd t_initial_condition)
    : OSNI::ODEb(t_polynomial_representation.getCoordinatesDimension(),
                 ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                 t_number_of_Chebyshev_points),
      m_B(t_polynomial_representation.m_B),
      m_strain_parameterisation(t_strain_parameterisation),
      m_internal_couples_integrator(t_internal_couples_integrator),
      m_internal_forces_integrator(t_internal_forces_integrator)
{
    this->setInitialConditions(t_initial_condition);
}

Eigen::VectorXd GeneralisedForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d C = m_internal_couples_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d N = m_internal_forces_integrator->getStateAtPoint(t_point);

    Eigen::Matrix<double, 6, 1> Lambda;
    Lambda << C,
              N;

    const auto Phi = m_strain_parameterisation->m_Phi_stack.at(t_point);
    return -Phi.transpose()*m_B.transpose()*Lambda;
}




IDMIntegrators::IDMIntegrators(const unsigned int t_number_of_Chebyshev_points,
                               const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                               std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                               std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_polynomial_representation(t_polynomial_representation),
      m_strain_parameterisation(t_strain_parameterisation),
      m_rod_properties(t_rod_properties)
{

}



}   //  namespace CROSP::IDMIntegrators
