 #include "CROSP/numerical_integrators/spectral_method/spectral_integrators.hpp"

namespace CROSP::numerical_integrators::spectral_method {


SpectralIntegrators::SpectralIntegrators(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                        std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                        std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
        : m_strain_parameterisation_stack(
              std::make_shared<::CROSP::strain_parameterisation_stack::StrainParameterisationStack>(t_strain_parameterisation->m_polynomial_representation,
                                                                                                    t_strain_parameterisation->m_number_of_Chebyshev_points)),
          m_Delta_strain_parameterisation_stack(
                    std::make_shared<::CROSP::strain_parameterisation_stack::StrainParameterisationStack>(t_strain_parameterisation->m_polynomial_representation,
                                                                                                          t_strain_parameterisation->m_number_of_Chebyshev_points,
                                                                                                          ::CROSP::strain_parameterisation_stack::zero_constrained_strain)),

          m_idm_integrators( std::make_shared<idm_integrators::IDMIntegrators>(m_strain_parameterisation_stack,
                                                                               t_rod_properties) ),
          m_tidm_integrators( /*std::make_shared<*/tidm_integrators::TIDMIntegrators/*>*/(m_strain_parameterisation_stack,
                                                                                  m_Delta_strain_parameterisation_stack,
                                                                                  m_idm_integrators,
                                                                                  t_rod_properties) )
    {}

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



    void SpectralIntegrators::updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                             const Eigen::VectorXd &t_Delta_dot_qe,
                                             const Eigen::VectorXd &t_Delta_ddot_qe)
    {
        m_Delta_strain_parameterisation_stack->updateStrainParameterisation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
    }

    void SpectralIntegrators::forwardTangentKinematics()
    {
//        //  Integrate Delta zeta
//        m_tidm_integrators->m_Delta_rotation->solveSystem();
//        m_tidm_integrators->m_Delta_position->solveSystem();

//        //  Integrate Delta eta
//        m_tidm_integrators->m_Delta_angular_velocity->solveSystem();
//        m_tidm_integrators->m_Delta_linear_velocity->solveSystem();

//        //  Integrate Delta dot eta
//        m_tidm_integrators->m_Delta_angular_acceleration->solveSystem();
//        m_tidm_integrators->m_Delta_linear_acceleration->solveSystem();
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



    void SpectralIntegrators::backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
    {

        Eigen::Vector3d Delta_couple_at_tip = t_Delta_Lambda_X1.block<3, 1>(0, 0);
        Eigen::Vector3d Delta_force_at_tip  = t_Delta_Lambda_X1.block<3, 1>(3, 0);


//        m_tidm_integrators->m_Delta_internal_forces->integrate(Delta_force_at_tip);
//        m_tidm_integrators->m_Delta_internal_couples->integrate(Delta_couple_at_tip);

//        m_tidm_integrators->m_Delta_generalised_forces->solveSystem();

    }




}   //  CROSP::numerical_integrators::spectral_method
