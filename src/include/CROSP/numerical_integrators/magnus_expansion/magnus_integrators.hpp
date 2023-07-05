#ifndef MAGNUS_INTEGRATORS_HPP
#define MAGNUS_INTEGRATORS_HPP

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"

//#include "magnus_idm_integrators/magnus_idm_integrators.hpp"

#include "CROSP/rod_properties/rod_properties.hpp"
#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"

#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"



namespace CROSP::numerical_integrators::magnus_expansion {


struct MagnusIntegrators {


    MagnusIntegrators(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                      std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                      std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);



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



    void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1);




    ::LieAlgebra::Vector6d getLambdaAtBase()const;



    LieAlgebra::Vector6d getDeltaLambdaAtBase()const;


    ::LieAlgebra::Vector6d getQaAtBase()const;



    LieAlgebra::Vector6d getDeltaQaAtBase()const;


    void updateIntegrationDomain(const double &t_rod_lenght);


    Eigen::MatrixXd getRodPositions()const;


    std::shared_ptr<::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_Delta_strain_parameterisation_stack;




    //  The set of integrators needed for the IDM



//    std::shared_ptr<magnus_idm_integrators::MagnusIDMIntegrators> m_inte;
    //magnus_idm_integrators::IDMIntegrators m_idm;
    //::CROSP::numerical_integrators::magnus_expansion::magnus_idm_integrators::IDMIntegrators m_idm;
//    ::CROSP::numerical_integrators::magnus_expansion::idm_integrators::MagnusIDMIntSPtr m_idm_integrators;
    //::CROSP::numerical_integrators::magnus_expansion::idm_integrators::MagnusIDMIntSPtr m_idm_integrators;

    //std::shared_ptr<CROSP::numerical_integrators::magnus_expansion::idm_integrators::IDMIntegrators> m_idm_integrators;

    //  The set of integrators needed for the TIDM
//    std::shared_ptr<tidm_integrators::TIDMIntegrators> m_tidm_integrators;



};





typedef std::shared_ptr<MagnusIntegrators> MagnusIntegratorsSPtr;






}   //  CROSP::numerical_integrators::magnus_expansion


#endif // MAGNUS_INTEGRATORS_HPP
