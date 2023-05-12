#ifndef SPECTRAL_INTEGRATORS_HPP
#define SPECTRAL_INTEGRATORS_HPP

#include "parameterisation_stack.hpp"
#include "idm_integrators.hpp"
#include "tidm_integrators.hpp"

namespace CROSP::spectral_integrators {


struct SpectralIntegrators {

    SpectralIntegrators(const strain_parameterisation::StrainParameterisation t_strain_parameterisation,
                        const strain_parameterisation::StrainParameterisation t_strain_parameterisation_Delta,
                        const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                        const rod_properties::RodProperties t_rod_properties,
                        const unsigned int t_number_of_Chebyshev_points)
        : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
          m_idm_integrators(
              std::make_shared<idm_integrators::IDMIntegrators>(t_strain_parameterisation,
                                                                t_polynomial_representation,
                                                                t_rod_properties,
                                                                t_number_of_Chebyshev_points)
              ),
          m_tidm_integrators(
              std::make_shared<tidm_integrators::TIDMIntegrators>(t_strain_parameterisation_Delta,
                                                                  t_polynomial_representation,
                                                                  m_idm_integrators,
                                                                  t_rod_properties,
                                                                  t_number_of_Chebyshev_points)
              )
    {}


    unsigned int m_number_of_Chebyshev_points;

    std::shared_ptr<idm_integrators::IDMIntegrators> m_idm_integrators;


    std::shared_ptr<tidm_integrators::TIDMIntegrators> m_tidm_integrators;


    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                const Eigen::VectorXd &t_dot_qe,
                                const Eigen::VectorXd &t_ddot_qe)
    {
        m_idm_integrators->m_parameterisation_stack->updateStacks(t_qe, t_dot_qe, t_ddot_qe);
    }


    void forwardKinematics()
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


    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
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


   void updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                    const Eigen::VectorXd &t_Delta_dot_qe,
                                    const Eigen::VectorXd &t_Delta_ddot_qe)
    {
        m_tidm_integrators->m_parameterisation_stack_Delta->updateStacks(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
    }

    void forwardTangentKinematics()
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


    void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
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






    ::LieAlgebra::Kinematics getKinematicsAtTip()const
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


    ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const
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



    void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1)
    {

        //  Map force and couple into local coordinates
        Eigen::Vector3d couple_at_tip = t_Lambda_X1.block<3, 1>(0, 0);
        Eigen::Vector3d force_at_tip  = t_Lambda_X1.block<3, 1>(3, 0);

        m_idm_integrators->m_internal_forces->integrate(force_at_tip);
        m_idm_integrators->m_internal_couples->integrate(couple_at_tip);


        m_idm_integrators->m_generalised_forces->solveSystem();

    }



    void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
    {

        Eigen::Vector3d Delta_couple_at_tip = t_Delta_Lambda_X1.block<3, 1>(0, 0);
        Eigen::Vector3d Delta_force_at_tip  = t_Delta_Lambda_X1.block<3, 1>(3, 0);


        m_tidm_integrators->m_Delta_internal_forces->integrate(Delta_force_at_tip);
        m_tidm_integrators->m_Delta_internal_couples->integrate(Delta_couple_at_tip);

        m_tidm_integrators->m_Delta_generalised_forces->solveSystem();
    }




    ::LieAlgebra::Vector6d getLambdaAtBase()const
    {
        ::LieAlgebra::Vector6d Lambda;
        Lambda <<   m_idm_integrators->m_internal_couples->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN),
                    m_idm_integrators->m_internal_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

        return Lambda;
    }



    LieAlgebra::Vector6d getDeltaLambdaAtBase()const
    {
        ::LieAlgebra::Vector6d Delta_Lambda;
        Delta_Lambda <<   m_tidm_integrators->m_Delta_internal_couples->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN),
                          m_tidm_integrators->m_Delta_internal_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

        return Delta_Lambda;
    }


    ::LieAlgebra::Vector6d getQaAtBase()const
    {
        return m_tidm_integrators->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);;
    }



    LieAlgebra::Vector6d getDeltaQaAtBase()const
    {
        return m_tidm_integrators->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);;
    }


    void updateIntegrationDomain(const double &t_rod_lenght)
    {
        m_idm_integrators->updateIntegrationDomain( t_rod_lenght );

        m_tidm_integrators->updateIntegrationDomain( t_rod_lenght );
    }


    Eigen::MatrixXd getRodPositions()const
    {
        return m_idm_integrators->m_position->getStackAsMatrix();
    }



};



}   //  namespace CROSP::spectral_integrator



#endif // SPECTRAL_INTEGRATORS_HPP
