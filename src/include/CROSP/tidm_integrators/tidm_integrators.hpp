/**
 * \file tidm_integrators.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod TIDM
 * \date 28-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#ifndef TIDM_INTEGRATORS_HPP
#define TIDM_INTEGRATORS_HPP

#include <memory>
#include <vector>

#include "OSNI/OSNI.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"


#include "CROSP/rod_properties/rod_properties.hpp"

#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/idm_integrators/idm_integrators.hpp"



namespace CROSP::tidm_integrators {









struct DeltaRotation : public OSNI::ODEAb {

    DeltaRotation(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                  std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                  const double &t_upper_integration_limit=1.0f,
                  const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Pointer to the perturbation of the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack_Delta;


    //  Storing pointers to direct access to the needed stacks
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
            m_parameterisation_stack->m_K_stack
        };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack {
        m_parameterisation_stack_Delta->m_K_stack
    };
};


struct DeltaPosition : public OSNI::ODEAb {

    DeltaPosition(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                  std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                  std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation,
                  const double &t_upper_integration_limit=1.0f,
                  const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Pointer to the perturbation of the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack_Delta;


    //  Storing pointers to direct access to the needed stacks
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_parameterisation_stack->m_Gamma_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Gamma_stack {
        m_parameterisation_stack_Delta->m_Gamma_stack
    };

    //  Pointers to the needed previous integrator(s)
    const std::shared_ptr<const OSNI::ODESolverInterface> m_delta_rotation;
};




struct DeltaAngularVelocity : public OSNI::ODEAb {

    DeltaAngularVelocity(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                         std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                         std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                         const double &t_upper_integration_limit=1.0f,
                         const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Pointer to the perturbation of the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack_Delta;


    //  Storing pointers to direct access to the needed stacks
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack {
        m_parameterisation_stack_Delta->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack {
        m_parameterisation_stack_Delta->m_dot_K_stack
    };


    //  Pointer to the IDM integrators for facilitate initialisation
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    //  Set of pointers to IDM integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator {
        m_idm_integrators->m_angular_velocity
    };

};


struct DeltaLinearVelocity : public OSNI::ODEAb {

    DeltaLinearVelocity(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                        std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                        std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                        std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                        const double &t_upper_integration_limit=1.0f,
                        const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;



    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Pointer to the perturbation of the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack_Delta;


    //  Storing pointers to direct access to the needed stacks
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack {
        m_parameterisation_stack_Delta->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_parameterisation_stack->m_Gamma_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Gamma_stack {
        m_parameterisation_stack_Delta->m_Gamma_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_Gamma_stack {
        m_parameterisation_stack_Delta->m_dot_Gamma_stack
    };


    //  Pointer to the IDM integrators for facilitate initialisation
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    //  Set of pointers to IDM integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator {
        m_idm_integrators->m_angular_velocity
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator {
        m_idm_integrators->m_linear_velocity
    };

    //  Set of pointers to the other Delta integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

};


struct DeltaAngularAccelerations : public OSNI::ODEAb {

    DeltaAngularAccelerations(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                              std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                              std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                              const double &t_upper_integration_limit=1.0f,
                              const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Pointer to the perturbation of the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack_Delta;


    //  Storing pointers to direct access to the needed stacks
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack {
        m_parameterisation_stack->m_dot_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack {
        m_parameterisation_stack_Delta->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack {
        m_parameterisation_stack_Delta->m_dot_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_ddot_K_stack {
        m_parameterisation_stack_Delta->m_ddot_K_stack
    };


    //  Pointer to the IDM integrators for facilitate initialisation
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    //  Set of pointers to IDM integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator {
        m_idm_integrators->m_angular_velocity
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration_integrator {
        m_idm_integrators->m_angular_acceleration
    };

    //  Set of pointers to the other Delta integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

};


struct DeltaLinearAccelerations : public OSNI::ODEAb {

    DeltaLinearAccelerations(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                             std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                             std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                             const double &t_upper_integration_limit=1.0f,
                             const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Pointer to the perturbation of the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack_Delta;


    //  Storing pointers to direct access to the needed stacks
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack {
        m_parameterisation_stack->m_dot_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack {
        m_parameterisation_stack_Delta->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack {
        m_parameterisation_stack_Delta->m_dot_K_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_parameterisation_stack->m_Gamma_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Gamma_stack {
        m_parameterisation_stack->m_dot_Gamma_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Gamma_stack {
        m_parameterisation_stack_Delta->m_Gamma_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_Gamma_stack {
        m_parameterisation_stack_Delta->m_dot_Gamma_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_ddot_Gamma_stack {
        m_parameterisation_stack_Delta->m_ddot_Gamma_stack
    };


    //  Pointer to the IDM integrators for facilitate initialisation
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    //  Set of pointers to IDM integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator {
        m_idm_integrators->m_angular_velocity
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator {
        m_idm_integrators->m_linear_velocity
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration_integrator {
        m_idm_integrators->m_angular_acceleration
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration_integrator{
        m_idm_integrators->m_linear_acceleration
    };

    //  Set of pointers to the other Delta integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_acceleration_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;

};




struct DeltaInternalForcesIntegrator : public OSNI::ODEAb {

    DeltaInternalForcesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                  std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                  std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                  const rod_properties::RodProperties t_rod_properties,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator,
                                  const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::Vector3d computeLocalExternalForces(unsigned int t_point) const;




    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Pointer to the perturbation of the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack_Delta;


    //  Storing pointers to direct access to the needed stacks
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack {
        m_parameterisation_stack_Delta->m_K_stack
    };


    //  Pointer to the IDM integrators for facilitate initialisation
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    //  Set of pointers to IDM integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator {
        m_idm_integrators->m_quaternion
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_position_integrator {
        m_idm_integrators->m_position
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator {
        m_idm_integrators->m_angular_velocity
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator {
        m_idm_integrators->m_linear_velocity
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator {
        m_idm_integrators->m_internal_forces
    };

    //  Set of pointers to the other Delta integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_rotation_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_acceleration_integrator;



    //  Pointer to the rod properties for the rod physics
    const rod_properties::RodProperties m_rod_properties;

    //  Store the linear block of the inertia matrix
    const Eigen::Matrix3d m_M_linear { m_rod_properties.getMLinear() };


};






struct DeltaInternalCouplesIntegrator : public OSNI::ODEAb {

    DeltaInternalCouplesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                   std::shared_ptr<const ParameterisationStack> t_parameterisation_stack_Delta,
                                   std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                   const rod_properties::RodProperties t_rod_properties,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                   const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;



    virtual Eigen::Vector3d computeLocalExternalCouples(unsigned int) const;



    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Pointer to the perturbation of the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack_Delta;


    //  Storing pointers to direct access to the needed stacks
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack {
        m_parameterisation_stack_Delta->m_K_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_parameterisation_stack->m_Gamma_stack
    };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Gamma_stack {
        m_parameterisation_stack_Delta->m_Gamma_stack
    };


    //  Pointer to the IDM integrators for facilitate initialisation
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    //  Set of pointers to IDM integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator {
        m_idm_integrators->m_quaternion
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_position_integrator {
        m_idm_integrators->m_position
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator {
        m_idm_integrators->m_angular_velocity
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator {
        m_idm_integrators->m_linear_velocity
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator {
        m_idm_integrators->m_internal_forces
    };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator {
        m_idm_integrators->m_internal_couples
    };


    //  Set of pointers to the other Delta integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_rotation_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_acceleration_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_forces_integrator;



    //  Pointer to the rod properties for the rod physics
    const rod_properties::RodProperties m_rod_properties;

    //  Store the two blocks of the inertia matrix separately
    const Eigen::Matrix3d m_M_angular { m_rod_properties.m_M.block<3,3>(0, 0) };
    const Eigen::Matrix3d m_M_linear { m_rod_properties.getMLinear() };

};



struct DeltaGeneralisedForcesIntegrator : public OSNI::ODEb {

    DeltaGeneralisedForcesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                     const double &t_upper_integration_limit=1.0f);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    //  Pointer to the strain parameterisation
    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    //  Set of pointers to the other Delta integrators
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_couples_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_forces_integrator;

};





struct TIDMIntegrators{

    TIDMIntegrators(const strain_parameterisation::StrainParameterisation t_strain_parameterisation_Delta,
                    const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                    const std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                    const rod_properties::RodProperties t_rod_properties,
                    const unsigned int t_number_of_Chebyshev_points);


    //  Instance of the rod properties
    rod_properties::RodProperties m_rod_properties;

    //  Instance of the variation of the strain parameterisation
    std::shared_ptr<ParameterisationStack> m_parameterisation_stack_Delta;

    //  Instance of the IDM integrators
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;






    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_rotation {
        std::make_shared<DeltaRotation>(m_idm_integrators->m_parameterisation_stack,
                                        m_parameterisation_stack_Delta,
                                        m_rod_properties.m_rod_dimensions.m_L)
    };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_position {
        std::make_shared<DeltaPosition>(m_idm_integrators->m_parameterisation_stack,
                                        m_parameterisation_stack_Delta,
                                        m_Delta_rotation,
                                        m_rod_properties.m_rod_dimensions.m_L)
    };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_angular_velocity {
        std::make_shared<DeltaAngularVelocity>(m_idm_integrators->m_parameterisation_stack,
                                               m_parameterisation_stack_Delta,
                                               m_idm_integrators,
                                               m_rod_properties.m_rod_dimensions.m_L)
    };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_linear_velocity {
        std::make_shared<DeltaLinearVelocity>(m_idm_integrators->m_parameterisation_stack,
                                              m_parameterisation_stack_Delta,
                                              m_idm_integrators,
                                              m_Delta_angular_velocity,
                                              m_rod_properties.m_rod_dimensions.m_L)
    };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_angular_acceleration {
        std::make_shared<DeltaAngularAccelerations>(m_idm_integrators->m_parameterisation_stack,
                                                    m_parameterisation_stack_Delta,
                                                    m_idm_integrators,
                                                    m_Delta_angular_velocity,
                                                    m_rod_properties.m_rod_dimensions.m_L)
    };




    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_linear_acceleration {
        std::make_shared<DeltaLinearAccelerations>(m_idm_integrators->m_parameterisation_stack,
                                                   m_parameterisation_stack_Delta,
                                                   m_idm_integrators,
                                                   m_Delta_angular_velocity,
                                                   m_Delta_angular_acceleration,
                                                   m_Delta_linear_velocity,
                                                   m_rod_properties.m_rod_dimensions.m_L)
    };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_forces {
        std::make_shared<DeltaInternalForcesIntegrator>(m_idm_integrators->m_parameterisation_stack,
                                                        m_parameterisation_stack_Delta,
                                                        m_idm_integrators,
                                                        m_rod_properties,
                                                        m_Delta_rotation,
                                                        m_Delta_angular_velocity,
                                                        m_Delta_linear_velocity,
                                                        m_Delta_linear_acceleration)
    };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_couples {
        std::make_shared<DeltaInternalCouplesIntegrator>(m_idm_integrators->m_parameterisation_stack,
                                                         m_parameterisation_stack_Delta,
                                                         m_idm_integrators,
                                                         m_rod_properties,
                                                         m_Delta_rotation,
                                                         m_Delta_angular_velocity,
                                                         m_Delta_angular_acceleration,
                                                         m_Delta_linear_velocity,
                                                         m_Delta_internal_forces)
    };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_generalised_forces {
        std::make_shared<DeltaGeneralisedForcesIntegrator>(m_idm_integrators->m_parameterisation_stack,
                                                           m_Delta_internal_couples,
                                                           m_Delta_internal_forces,
                                                           m_rod_properties.m_rod_dimensions.m_L)
    };




    /*!
     * \brief updateIntegrationDomain this function is used to set the integration upper limit at every integrator
     * \param t_upper_integration_limit is the upper integration limit used to scale the normalised domain [0, 1]
     */
    void updateIntegrationDomain(const double &t_upper_integration_limit=1.0f);
};

}   //  namespace CROSP

#endif // TIDM_INTEGRATORS_HPP
