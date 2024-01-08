#ifndef SPECTRAL_INTEGRATORS_HPP
#define SPECTRAL_INTEGRATORS_HPP

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"

#include "idm_integrators/idm_integrators.hpp"
#include "tidm_integrators/tidm_integrators.hpp"
#include "internal_actuation/internal_actuation_integrator.hpp"
#include "integration_polynomial_base/integration_polynomial_base.hpp"



namespace CROSP::numerical_integrators::spectral_method {


struct SpectralIntegrators{


    SpectralIntegrators(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                        const unsigned int t_number_of_Chebyshev_points,
                        std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                        strain_parameterisation_stack::StrainFunction t_Xi_c=::CROSP::strain_parameterisation_stack::default_constrained_strain);

//    SpectralIntegrators(SpectralIntegrators &t_other);


    void addInternalActuation(::ATORS::distributed_actuation::DistributedActuationUptr t_distributed_actuation);


    std::string printIntegratorProperties()const;


    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                        const Eigen::VectorXd &t_dot_qe,
                                        const Eigen::VectorXd &t_ddot_qe);


    void forwardKinematics();


    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                   const Eigen::Vector3d &t_initial_position,
                                   const Eigen::Vector3d &t_initial_angular_velocity,
                                   const Eigen::Vector3d &t_initial_linear_velocity,
                                   const Eigen::Vector3d &t_initial_angular_acceleration,
                                   const Eigen::Vector3d &t_initial_linear_acceleration);


   void updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                            const Eigen::VectorXd &t_Delta_dot_qe,
                                            const Eigen::VectorXd &t_Delta_ddot_qe);

    void forwardTangentKinematics();


    void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                          const Eigen::Vector3d &t_initial_Delta_position,
                                          const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                          const Eigen::Vector3d &t_initial_Delta_linear_acceleration);


    ::LieAlgebra::Kinematics getKinematicsAtTip()const;

    ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const;



    void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1);


    void updateInternalActuation(const double &t_time);

    Eigen::VectorXd getQad() const;



    void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1);




    ::LieAlgebra::Vector6d getLambdaAtBase()const;



    LieAlgebra::Vector6d getDeltaLambdaAtBase()const;


    Eigen::VectorXd getQaAtBase()const;



    Eigen::VectorXd getDeltaQaAtBase()const;


    void updateIntegrationDomain(const double &t_rod_lenght);


    Eigen::MatrixXd getRodPositions()const;


    FullODEStatesObservations getIDMStatesObservations();

    FullODEStatesObservations getTIDMStatesObservations();


    Eigen::MatrixXd integratePhiTPhi()const;


    std::shared_ptr<::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_Delta_strain_parameterisation_stack;




    //  The set of integrators needed for the IDM
    std::shared_ptr<idm_integrators::IDMIntegrators> m_idm_integrators;

    //  The set of integrators needed for the TIDM
    std::shared_ptr<tidm_integrators::TIDMIntegrators> m_tidm_integrators;


    internal_actuation_integrator::InternalActuationIntegratorUptr m_internal_actuation_integrator {
        std::make_unique<internal_actuation_integrator::NullInternalActuationIntegrator>(m_strain_parameterisation_stack)
    };


    [[deprecated("This member should be deleted as it is not safe to do the integration at once but is should be done with Hcal")]]
    std::unique_ptr<integration_polynomial_base::IntegratorPolynomialBase> m_integrator_polynomial_base {
        std::make_unique<integration_polynomial_base::IntegratorPolynomialBase>(m_strain_parameterisation_stack)
    };



};




typedef std::shared_ptr<SpectralIntegrators> IntegratorsSPtr;






}   //  CROSP::numerical_integrators::spectral_method

#endif // SPECTRAL_INTEGRATORS_HPP
