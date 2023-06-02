/**
 * \file idm_odes.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the odes for the Cosserat rod IDM
 * \date 1-06-2023
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#ifndef IDM_ODES_HPP
#define IDM_ODES_HPP

#include <memory>

#include "CROSP/rod_properties/rod_properties.hpp"
#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include "OSNI/OSNI.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"

#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"


/// \brief CROSP::idm_integrators is the namespace containing the integrators for the IDM
namespace CROSP::numerical_integrators::spectral_method::idm_odes {







/*!
 * \brief The QuaternionIntegrator struct implements the spectral numerical integration for the quaternion
 *
 * This class implemnts the integration for the quaternions forllowing the formula \f$ Q' = \frac{1}{2}  A_{(K)} Q
 */
struct QuaternionODE {

    QuaternionODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack);

    Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point);

    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_strain_parameterisation_stack->m_K_stack
    };

    Eigen::Matrix4d m_A_at_chebychev_point;

};



struct PositionODE {

    PositionODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator);


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_strain_parameterisation_stack->m_Gamma_stack
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    Eigen::Quaterniond m_quaternion_at_point;

};




struct AngularVelocityODE {

    AngularVelocityODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack {
        m_strain_parameterisation_stack->m_dot_K_stack
    };


};


struct LinearVelocityODE {

    LinearVelocityODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                      std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point);


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point);


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_Gamma_stack {
        m_strain_parameterisation_stack->m_hat_Gamma_stack
    };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Gamma_stack {
        m_strain_parameterisation_stack->m_dot_Gamma_stack
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;



};




struct AngularAccelerationODE {

    AngularAccelerationODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                           std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_K_stack {
        m_strain_parameterisation_stack->m_hat_dot_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_K_stack {
        m_strain_parameterisation_stack->m_ddot_K_stack
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;


};







struct LinearAccelerationODE {

    LinearAccelerationODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                          std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                          std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                          std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator);



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_K_stack {
        m_strain_parameterisation_stack->m_hat_dot_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_Gamma_stack {
        m_strain_parameterisation_stack->m_hat_Gamma_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_Gamma_stack{
        m_strain_parameterisation_stack->m_hat_dot_Gamma_stack
    };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_Gamma_stack{
        m_strain_parameterisation_stack->m_ddot_Gamma_stack
    };

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;



};




struct InternalForcesODE {

    InternalForcesODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                      std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                      std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                      std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                      std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                      std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                      std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator);



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computeDistributedForce(const unsigned int t_point) const;




    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position;

};


struct InternalCouplesODE {

    InternalCouplesODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                       std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator);




    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computeDistributedCouple(const unsigned int) const;


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_strain_parameterisation_stack->m_Gamma_stack
    };

    std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position;

    std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces;

};


struct GeneralisedForcesODE {

    GeneralisedForcesODE(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                         std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
                         std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    unsigned int m_coordinates_dimension;

    std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator;
    std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;

};




}   //  namespace CROSP::numerical_integrators::spectral_method::idm_odes





#endif // IDM_ODES_HPP
