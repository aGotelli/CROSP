#ifndef SPECTRAL_INTEGRATORS_HPP
#define SPECTRAL_INTEGRATORS_HPP

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"

#include "idm_integrators/idm_integrators.hpp"
#include "tidm_integrators/tidm_integrators.hpp"

namespace CROSP::numerical_integrators::spectral_method {


struct SpectralIntegrators : public CosseratRodIntegrators {


    SpectralIntegrators(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                        const unsigned int t_number_of_Chebyshev_points,
                        std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);



    virtual std::string printIntegratorProperties()const final;


    virtual void updateParameterisation(const Eigen::VectorXd &t_qe,
                                        const Eigen::VectorXd &t_dot_qe,
                                        const Eigen::VectorXd &t_ddot_qe) final;


    virtual void forwardKinematics() final;


    virtual void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                   const Eigen::Vector3d &t_initial_position,
                                   const Eigen::Vector3d &t_initial_angular_velocity,
                                   const Eigen::Vector3d &t_initial_linear_velocity,
                                   const Eigen::Vector3d &t_initial_angular_acceleration,
                                   const Eigen::Vector3d &t_initial_linear_acceleration) final;


   virtual void updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                            const Eigen::VectorXd &t_Delta_dot_qe,
                                            const Eigen::VectorXd &t_Delta_ddot_qe) final;

    virtual void forwardTangentKinematics() final;


    virtual void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                          const Eigen::Vector3d &t_initial_Delta_position,
                                          const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                          const Eigen::Vector3d &t_initial_Delta_linear_acceleration) final;


    virtual ::LieAlgebra::Kinematics getKinematicsAtTip()const final;

    virtual ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const final;



    virtual void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1) final;



    virtual void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1) final;




    virtual ::LieAlgebra::Vector6d getLambdaAtBase()const final;



    virtual LieAlgebra::Vector6d getDeltaLambdaAtBase()const final;


    virtual ::LieAlgebra::Vector6d getQaAtBase()const final;



    virtual LieAlgebra::Vector6d getDeltaQaAtBase()const final;


    virtual void updateIntegrationDomain(const double &t_rod_lenght) final;


    virtual Eigen::MatrixXd getRodPositions()const final;


    std::shared_ptr<::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_Delta_strain_parameterisation_stack;




    //  The set of integrators needed for the IDM
    std::shared_ptr<idm_integrators::IDMIntegrators> m_idm_integrators;

    //  The set of integrators needed for the TIDM
    std::shared_ptr<tidm_integrators::TIDMIntegrators> m_tidm_integrators;



};


//SpectralIntegrators::SpectralIntegrators(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
//                                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
//                                         std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
//    :   m_idm_integrators( std::make_shared<idm_integrators::IDMIntegrators>(t_strain_parameterisation,
//                                                                             t_rod_properties) ),
//        m_tidm_integrators( std::make_shared<tidm_integrators::TIDMIntegrators>(t_strain_parameterisation,
//                                                                                t_strain_parameterisation_Delta,
//                                                                                m_idm_integrators,
//                                                                                t_rod_properties) )
//{}



typedef std::shared_ptr<SpectralIntegrators> IntegratorsSPtr;






}   //  CROSP::numerical_integrators::spectral_method

#endif // SPECTRAL_INTEGRATORS_HPP
