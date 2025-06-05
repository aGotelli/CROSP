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


#include "CROSP/numerical_integrators/spectral_method/idm_integrators/idm_integrators.hpp"



namespace CROSP::numerical_integrators::spectral_method::idm_integrators {






QuaternionIntegrator::QuaternionIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                           const double &t_upper_integration_limit,
                                           const Eigen::Vector4d &t_initial_condition)
    : OSNI::ODEA(4,
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_strain_parameterisation_stack->m_number_of_points,
                 t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack)
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







AngularVelocityIntegrator::AngularVelocityIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                     const double &t_upper_integration_limit,
                                                     const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack)
{
    this->setInitialConditions(t_initial_condition);
}


Eigen::MatrixXd AngularVelocityIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}


Eigen::VectorXd AngularVelocityIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_dot_K_stack->at(t_point);
}





LinearVelocityIntegrator::LinearVelocityIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                   const double &t_upper_integration_limit,
                                                   const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_angular_velocity(t_angular_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}


Eigen::MatrixXd LinearVelocityIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}


Eigen::VectorXd LinearVelocityIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d dot_Gamma = m_dot_Gamma_stack->at(t_point);
    Eigen::Vector3d Omega = m_angular_velocity->getStateAtPoint(t_point);

    return dot_Gamma - m_hat_Gamma_stack->at(t_point)*Omega;
}





AngularAccelerationIntegrator::AngularAccelerationIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                             const double &t_upper_integration_limit,
                                                             const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_angular_velocity(t_angular_velocity_integrator)
{
    this->setInitialConditions(t_initial_condition);
}



Eigen::MatrixXd AngularAccelerationIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}


Eigen::VectorXd AngularAccelerationIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_ddot_K_stack->at(t_point)
            - m_hat_dot_K_stack->at(t_point) * m_angular_velocity->getStateAtPoint(t_point);
}









LinearAccelerationIntegrator::LinearAccelerationIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                           std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                                           const double &t_upper_integration_limit,
                                                           const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_angular_velocity(t_angular_velocity_integrator),
      m_linear_velocity(t_linear_velocity_integrator),
      m_angular_acceleration(t_angular_acceleration_integrator)
{
    this->setInitialConditions(t_initial_condition);
}



Eigen::MatrixXd LinearAccelerationIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    return - m_hat_K_stack->at(t_point);
}


Eigen::VectorXd LinearAccelerationIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    return m_ddot_Gamma_stack->at(t_point)
            - m_hat_Gamma_stack->at(t_point) * m_angular_acceleration->getStateAtPoint(t_point)
            - m_hat_dot_Gamma_stack->at(t_point) * m_angular_velocity->getStateAtPoint(t_point)
            - m_hat_dot_K_stack->at(t_point) * m_linear_velocity->getStateAtPoint(t_point);
}






InternalForcesIntegrator::InternalForcesIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                   std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                                                   std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                                                   const double &t_upper_integration_limit,
                                                   const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_rod_properties(t_rod_properties),
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
    return  m_hat_K_stack->at(t_point).transpose();
}


Eigen::VectorXd InternalForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    const auto M_linear = m_rod_properties->getMLinear();


    Eigen::Vector3d inertial_acceleration = M_linear*m_linear_acceleration->getStateAtPoint(t_point);
    Eigen::Vector3d inertial_velocities = ::LieAlgebra::skew( m_angular_velocity->getStateAtPoint(t_point) ).transpose() * M_linear * m_linear_velocity->getStateAtPoint(t_point);
    Eigen::Vector3d distributed_forces = computeDistributedForce(t_point);


    // std::cout << "At point : " << t_point << "\n" <<
    //     "\t dot_V : \n" << m_linear_acceleration->getStateAtPoint(t_point) << "\n" << std::endl;
    //     // "\t inertial_velocities : \n" << inertial_velocities << "\n" <<
    //     // "\t distributed_forces : \n" << distributed_forces << std::endl;

    Eigen::Vector3d b = inertial_acceleration
                        - inertial_velocities
                        - distributed_forces
                        ;


    return b;
}


Eigen::VectorXd InternalForcesIntegrator::computeDistributedForce(const unsigned int t_point) const
{
//    if(t_point == 0 or t_point == m_strain_parameterisation_stack->m_number_of_points-1)
//        return Eigen::Vector3d::Zero();

    Eigen::Vector4d q(m_quaternion->getStateAtPoint(t_point));

    Eigen::Matrix3d R = Eigen::Quaterniond(q(0),
                                           q(1),
                                           q(2),
                                           q(3)).toRotationMatrix();


    const Eigen::Vector3d distributed_weight_force = m_rod_properties->distributedGravitationalForce();
    Eigen::Vector3d N_bar = R.transpose()*distributed_weight_force;

    return N_bar;
}




InternalCouplesIntegrator::InternalCouplesIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                     std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                                                     std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                                                     const double &t_upper_integration_limit,
                                                     const Eigen::Vector3d &t_initial_condition)
    : OSNI::ODEAb(3,
                  ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                  t_strain_parameterisation_stack->m_number_of_points,
                  t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_rod_properties( t_rod_properties ),
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
    return  m_hat_K_stack->at(t_point).transpose();
}


Eigen::VectorXd InternalCouplesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{

    const auto M_linear = m_rod_properties->getMLinear();
    const auto M_angular = m_rod_properties->getMAngular();

    return ::LieAlgebra::skew( m_Gamma_stack->at(t_point) ).transpose()*m_internal_forces->getStateAtPoint(t_point)
            + M_angular*m_angular_acceleration->getStateAtPoint(t_point)
            - ::LieAlgebra::skew(m_angular_velocity->getStateAtPoint(t_point)).transpose() * M_angular * m_angular_velocity->getStateAtPoint(t_point)
            - ::LieAlgebra::skew( m_linear_velocity->getStateAtPoint(t_point) ).transpose()* M_linear * m_linear_velocity->getStateAtPoint(t_point)
            + computeDistributedCouple(t_point);
}


Eigen::VectorXd InternalCouplesIntegrator::computeDistributedCouple(const unsigned int) const
{
    return Eigen::Vector3d::Zero();
}




GeneralisedForcesIntegrator::GeneralisedForcesIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                         std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
                                                         std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                                                         const double &t_upper_integration_limit)
    : OSNI::ODEb(t_strain_parameterisation_stack->m_B_Phi_stack.cols(),
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_strain_parameterisation_stack->m_number_of_points,
                 t_upper_integration_limit),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_coordinates_dimension(t_strain_parameterisation_stack->m_B_Phi_stack.cols()),
      m_internal_couples_integrator(t_internal_couples_integrator),
      m_internal_forces_integrator(t_internal_forces_integrator)
{
    this->setInitialConditions( Eigen::VectorXd::Zero(m_coordinates_dimension) );
}

Eigen::VectorXd GeneralisedForcesIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    Eigen::Vector3d C = m_internal_couples_integrator->getStateAtPoint(t_point);
    Eigen::Vector3d N = m_internal_forces_integrator->getStateAtPoint(t_point);

    Eigen::Matrix<double, 6, 1> Lambda;
    Lambda << C,
              N;

    const Eigen::VectorXd b_at_point = m_strain_parameterisation_stack->m_B_Phi_stack.block(6*t_point, 0, 6, m_coordinates_dimension).transpose()*Lambda;
    return b_at_point;
}




IDMIntegrators::IDMIntegrators(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                               std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
    : m_strain_parameterisation_stack( t_strain_parameterisation_stack ),
      m_rod_properties(t_rod_properties)
{}

void IDMIntegrators::updateIntegrationDomain(const double &t_upper_integration_limit)
{
    m_quaternion->setUpperIntegrationDomain(t_upper_integration_limit);
    m_position->setUpperIntegrationDomain(t_upper_integration_limit);

    m_angular_velocity->setUpperIntegrationDomain(t_upper_integration_limit);
    m_linear_velocity->setUpperIntegrationDomain(t_upper_integration_limit);

    m_angular_acceleration->setUpperIntegrationDomain(t_upper_integration_limit);
    m_linear_acceleration->setUpperIntegrationDomain(t_upper_integration_limit);

    m_internal_forces->setUpperIntegrationDomain(t_upper_integration_limit);
    m_internal_couples->setUpperIntegrationDomain(t_upper_integration_limit);

    m_generalised_forces->setUpperIntegrationDomain(t_upper_integration_limit);
}



}   //  namespace CROSP::numerical_integrators::spectral_method::idm_integrators
