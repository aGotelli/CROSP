#ifndef MAGNUS_INTEGRATORS_HPP
#define MAGNUS_INTEGRATORS_HPP

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"

#include "ATORS/ATORS.hpp"

#include "CROSP/rod_properties/rod_properties.hpp"
#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"

#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"


#include "magnus_idm_integrators/magnus_idm_integrators.hpp"

namespace CROSP::numerical_integrators::magnus_expansion {


struct MagnusIntegrators {


    MagnusIntegrators(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                      std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                      std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);


    void addInternalActuation(::ATORS::distributed_actuation::DistributedActuationUptr t_distributed_actuation)
    {}

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


    Eigen::VectorXd getQaAtBase()const;



    Eigen::VectorXd getDeltaQaAtBase()const;


    void updateIntegrationDomain(const double &t_rod_lenght);


    Eigen::MatrixXd getRodPositions()const;


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_Delta_strain_parameterisation_stack;


    std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;


    std::shared_ptr<magnus_expansion::magnus_idm_integrators::MagnusIDMIntegrators> m_idm_integrators {
        std::make_shared<magnus_expansion::magnus_idm_integrators::MagnusIDMIntegrators>(m_strain_parameterisation_stack,
                                                                                          m_rod_properties)
    };



};





typedef std::shared_ptr<MagnusIntegrators> MagnusIntegratorsSPtr;






}   //  CROSP::numerical_integrators::magnus_expansion


#endif // MAGNUS_INTEGRATORS_HPP
