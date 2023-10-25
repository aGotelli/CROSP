 #include "CROSP/numerical_integrators/spectral_method/spectral_integrators.hpp"

namespace CROSP::numerical_integrators::spectral_method {


SpectralIntegrators::SpectralIntegrators(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                         const unsigned int t_number_of_Chebyshev_points,
                                         std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                         ::CROSP::strain_parameterisation_stack::StrainFunction t_Xi_c)
    : m_strain_parameterisation_stack(
          std::make_shared<::CROSP::strain_parameterisation_stack::StrainParameterisationStack>(t_polynomial_representation,
                                                                                                t_number_of_Chebyshev_points,
                                                                                                t_Xi_c)),
      m_Delta_strain_parameterisation_stack(
          std::make_shared<::CROSP::strain_parameterisation_stack::StrainParameterisationStack>(t_polynomial_representation,
                                                                                                t_number_of_Chebyshev_points,
                                                                                                ::CROSP::strain_parameterisation_stack::zero_constrained_strain)),

      m_idm_integrators( std::make_shared<idm_integrators::IDMIntegrators>(m_strain_parameterisation_stack,
                                                                           t_rod_properties) ),
      m_tidm_integrators( std::make_shared<tidm_integrators::TIDMIntegrators>(m_strain_parameterisation_stack,
                                                                              m_Delta_strain_parameterisation_stack,
                                                                              m_idm_integrators,
                                                                              t_rod_properties) )
{}


void SpectralIntegrators::addInternalActuation(::ATORS::distributed_actuation::DistributedActuationUptr t_distributed_actuation)
{
    m_internal_actuation_integrator =
            std::make_unique<internal_actuation_integrator::InternalActuationIntegrator>(m_strain_parameterisation_stack,
                                                                                         std::move(t_distributed_actuation),
                                                                                         m_idm_integrators->m_rod_properties->m_rod_dimensions.m_L);
}


std::string SpectralIntegrators::printIntegratorProperties()const
{
    std::stringstream integrator_properties;

    integrator_properties << "      Spectral integrator\n"
                             "          Number of Chebyshev points : " << m_strain_parameterisation_stack->m_number_of_points << "\n"
                             "          Integration domain : [0, 1]\n";

    return integrator_properties.str();

}

void SpectralIntegrators::updateParameterisation(const Eigen::VectorXd &t_qe,
                                                 const Eigen::VectorXd &t_dot_qe,
                                                 const Eigen::VectorXd &t_ddot_qe)
{
    m_strain_parameterisation_stack->updateStrainParameterisation(t_qe, t_dot_qe, t_ddot_qe);
}


void SpectralIntegrators::forwardKinematics()
{
    //  Integrate Quaternions
    m_idm_integrators->m_quaternion->solveSystem();

    //  Integrate Positions
    m_idm_integrators->m_position->solveSystem();

    //  Integrate twist
    m_idm_integrators->m_angular_velocity->solveSystem();
    m_idm_integrators->m_linear_velocity->solveSystem();

    //  Integrate accelerations
    m_idm_integrators->m_angular_acceleration->solveSystem();
    m_idm_integrators->m_linear_acceleration->solveSystem();
}


void SpectralIntegrators::forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                            const Eigen::Vector3d &t_initial_position,
                                            const Eigen::Vector3d &t_initial_angular_velocity,
                                            const Eigen::Vector3d &t_initial_linear_velocity,
                                            const Eigen::Vector3d &t_initial_angular_acceleration,
                                            const Eigen::Vector3d &t_initial_linear_acceleration)
{
    //  Integrate Quaternions
    m_idm_integrators->m_quaternion->integrate(t_initial_quaternion);

    //  Integrate Positions
    m_idm_integrators->m_position->integrate(t_initial_position);

    //  Integrate twist
    m_idm_integrators->m_angular_velocity->integrate(t_initial_angular_velocity);
    m_idm_integrators->m_linear_velocity->integrate(t_initial_linear_velocity);

    //  Integrate accelerations
    m_idm_integrators->m_angular_acceleration->integrate(t_initial_angular_acceleration);
    m_idm_integrators->m_linear_acceleration->integrate(t_initial_linear_acceleration);
}



void SpectralIntegrators::updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                                      const Eigen::VectorXd &t_Delta_dot_qe,
                                                      const Eigen::VectorXd &t_Delta_ddot_qe)
{
    m_Delta_strain_parameterisation_stack->updateStrainParameterisation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
}

void SpectralIntegrators::forwardTangentKinematics()
{
    //  Integrate Delta zeta
    m_tidm_integrators->m_Delta_rotation->solveSystem();
    m_tidm_integrators->m_Delta_position->solveSystem();

    //  Integrate Delta eta
    m_tidm_integrators->m_Delta_angular_velocity->solveSystem();
    m_tidm_integrators->m_Delta_linear_velocity->solveSystem();

    //  Integrate Delta dot eta
    m_tidm_integrators->m_Delta_angular_acceleration->solveSystem();
    m_tidm_integrators->m_Delta_linear_acceleration->solveSystem();
}



void SpectralIntegrators::forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                                   const Eigen::Vector3d &t_initial_Delta_position,
                                                   const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                                   const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                                   const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                                   const Eigen::Vector3d &t_initial_Delta_linear_acceleration)
{
    //  Integrate Delta zeta
    m_tidm_integrators->m_Delta_rotation->integrate( t_initial_Delta_orientation );
    m_tidm_integrators->m_Delta_position->integrate( t_initial_Delta_position );

    //  Integrate Delta eta
    m_tidm_integrators->m_Delta_angular_velocity->integrate( t_initial_Delta_angular_velocity );
    m_tidm_integrators->m_Delta_linear_velocity->integrate( t_initial_Delta_linear_velocity );

    //  Integrate Delta dot eta
    m_tidm_integrators->m_Delta_angular_acceleration->integrate( t_initial_Delta_angular_acceleration );
    m_tidm_integrators->m_Delta_linear_acceleration->integrate( t_initial_Delta_linear_acceleration );
}






::LieAlgebra::Kinematics SpectralIntegrators::getKinematicsAtTip()const
{
    ::LieAlgebra::Kinematics rod_tip_kinematics;

    rod_tip_kinematics.m_pose =
            ::LieAlgebra::SE3Pose( m_idm_integrators->m_quaternion->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
                                   m_idm_integrators->m_position->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END) );


    rod_tip_kinematics.m_twist <<
                                  m_idm_integrators->m_angular_velocity->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
            m_idm_integrators->m_linear_velocity->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

    rod_tip_kinematics.m_accelerations <<
                                          m_idm_integrators->m_angular_acceleration->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
            m_idm_integrators->m_linear_acceleration->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

    return rod_tip_kinematics;
}


::LieAlgebra::TangentKinematics SpectralIntegrators::getTangentKinematicsAtTip()const
{
    ::LieAlgebra::TangentKinematics rod_tip_tangent_kinematics;

    rod_tip_tangent_kinematics.m_Delta_zeta <<
                                               m_tidm_integrators->m_Delta_rotation->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
            m_tidm_integrators->m_Delta_position->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

    rod_tip_tangent_kinematics.m_Delta_twist <<
                                                m_tidm_integrators->m_Delta_angular_velocity->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
            m_tidm_integrators->m_Delta_linear_velocity->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

    rod_tip_tangent_kinematics.m_Delta_acceleration <<
                                                       m_tidm_integrators->m_Delta_angular_acceleration->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
            m_tidm_integrators->m_Delta_linear_acceleration->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

    return rod_tip_tangent_kinematics;
}






void SpectralIntegrators::backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1)
{

    //  Map force and couple into local coordinates
    Eigen::Vector3d couple_at_tip = t_Lambda_X1.block<3, 1>(0, 0);
    Eigen::Vector3d force_at_tip  = t_Lambda_X1.block<3, 1>(3, 0);

    m_idm_integrators->m_internal_forces->integrate(force_at_tip);
    m_idm_integrators->m_internal_couples->integrate(couple_at_tip);


    m_idm_integrators->m_generalised_forces->solveSystem();

}




void SpectralIntegrators::updateInternalActuation(const double &t_time)
{
    m_internal_actuation_integrator->updateActuation(t_time);


}


Eigen::VectorXd SpectralIntegrators::getQad()const
{
    const Eigen::VectorXd Q_ad = m_internal_actuation_integrator->getActuation();

    return Q_ad;
}



void SpectralIntegrators::backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
{

    Eigen::Vector3d Delta_couple_at_tip = t_Delta_Lambda_X1.block<3, 1>(0, 0);
    Eigen::Vector3d Delta_force_at_tip  = t_Delta_Lambda_X1.block<3, 1>(3, 0);


    m_tidm_integrators->m_Delta_internal_forces->integrate(Delta_force_at_tip);
    m_tidm_integrators->m_Delta_internal_couples->integrate(Delta_couple_at_tip);

    m_tidm_integrators->m_Delta_generalised_forces->solveSystem();

}



::LieAlgebra::Vector6d SpectralIntegrators::getLambdaAtBase()const
{
    ::LieAlgebra::Vector6d Lambda;
    Lambda <<   m_idm_integrators->m_internal_couples->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN),
            m_idm_integrators->m_internal_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

    return Lambda;
}



LieAlgebra::Vector6d SpectralIntegrators::getDeltaLambdaAtBase()const
{
    ::LieAlgebra::Vector6d Delta_Lambda;
    Delta_Lambda <<   m_tidm_integrators->m_Delta_internal_couples->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN),
            m_tidm_integrators->m_Delta_internal_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

    return Delta_Lambda;
}


::LieAlgebra::Vector6d SpectralIntegrators::getQaAtBase()const
{
    return m_tidm_integrators->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);;
}



LieAlgebra::Vector6d SpectralIntegrators::getDeltaQaAtBase()const
{
    return m_tidm_integrators->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);;
}


void SpectralIntegrators::updateIntegrationDomain(const double &t_rod_lenght)
{
    m_idm_integrators->updateIntegrationDomain( t_rod_lenght );

    m_tidm_integrators->updateIntegrationDomain( t_rod_lenght );
}


Eigen::MatrixXd SpectralIntegrators::getRodPositions()const
{
    return m_idm_integrators->m_position->getStackAsMatrix();
}


FullODEStatesObservations SpectralIntegrators::getIDMStatesObservations()
{
    FullODEStatesObservations ODE_states_observations;


    ODE_states_observations.orientation_stack = m_idm_integrators->m_quaternion->getStackAsMatrix();
    ODE_states_observations.r_stack = m_idm_integrators->m_position->getStackAsMatrix();

    ODE_states_observations.Omega_stack = m_idm_integrators->m_angular_velocity->getStackAsMatrix();
    ODE_states_observations.V_stack = m_idm_integrators->m_linear_velocity->getStackAsMatrix();

    ODE_states_observations.dot_Omega_stack = m_idm_integrators->m_angular_acceleration->getStackAsMatrix();
    ODE_states_observations.dot_V_stack = m_idm_integrators->m_linear_acceleration->getStackAsMatrix();

    ODE_states_observations.C_stack = m_idm_integrators->m_internal_couples->getStackAsMatrix();
    ODE_states_observations.N_stack = m_idm_integrators->m_internal_forces->getStackAsMatrix();

    ODE_states_observations.Qa_stack = m_idm_integrators->m_generalised_forces->getStackAsMatrix();
    ODE_states_observations.Qad_stack = m_internal_actuation_integrator->getStackAsMatrix();

    return ODE_states_observations;

}


FullODEStatesObservations SpectralIntegrators::getTIDMStatesObservations()
{
    FullODEStatesObservations ODE_states_observations;


    ODE_states_observations.orientation_stack = m_tidm_integrators->m_Delta_rotation->getStackAsMatrix();
    ODE_states_observations.r_stack = m_tidm_integrators->m_Delta_position->getStackAsMatrix();

    ODE_states_observations.Omega_stack = m_tidm_integrators->m_Delta_angular_velocity->getStackAsMatrix();
    ODE_states_observations.V_stack = m_tidm_integrators->m_Delta_linear_velocity->getStackAsMatrix();

    ODE_states_observations.dot_Omega_stack = m_tidm_integrators->m_Delta_angular_acceleration->getStackAsMatrix();
    ODE_states_observations.dot_V_stack = m_tidm_integrators->m_Delta_linear_acceleration->getStackAsMatrix();

    ODE_states_observations.C_stack = m_tidm_integrators->m_Delta_internal_couples->getStackAsMatrix();
    ODE_states_observations.N_stack = m_tidm_integrators->m_Delta_internal_forces->getStackAsMatrix();

    ODE_states_observations.Qa_stack = m_tidm_integrators->m_Delta_generalised_forces->getStackAsMatrix();

    return ODE_states_observations;

}

Eigen::MatrixXd SpectralIntegrators::integratePhiTPhi()const
{
    m_integrator_polynomial_base->solveSystem();

    Eigen::VectorXd PhiTPhi_state = m_integrator_polynomial_base->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

    unsigned int size = sqrt(PhiTPhi_state.size());

    Eigen::MatrixXd PhiTPhi = Eigen::Map<Eigen::MatrixXd>(PhiTPhi_state.data(), size, size);

    return PhiTPhi;
}



}   //  CROSP::numerical_integrators::spectral_method
