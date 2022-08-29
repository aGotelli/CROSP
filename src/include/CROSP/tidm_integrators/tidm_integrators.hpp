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

#include "CROSP/rod_properties/rod_properties.hpp"

#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/idm_integrators/idm_integrators.hpp"

namespace CROSP::tidm_integrators {



struct DeltaRotation : public OSNI::ODEAb {
    DeltaRotation(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
};


struct DeltaPosition : public OSNI::ODEAb {

    DeltaPosition(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                  std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_delta_rotation;
};




struct DeltaAngularVelocity : public OSNI::ODEAb {
    DeltaAngularVelocity(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                         std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;

};


struct DeltaLinearVelocity : public OSNI::ODEAb {
    DeltaLinearVelocity(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                        std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                        std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                        std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;



    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;

};


struct DeltaAngularAccelerations : public OSNI::ODEAb {

    DeltaAngularAccelerations(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                              std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                              std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_ddot_K_stack;


    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

};


struct DeltaLinearAccelerations : public OSNI::ODEAb {

    DeltaLinearAccelerations(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                             std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_ddot_Lambda_stack;


    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_acceleration_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;

};




struct DeltaInternalForcesIntegrator : public OSNI::ODEAb {

    DeltaInternalForcesIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                  std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                  std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::Vector3d computeLocalExternalForces(unsigned int t_point) const;



    const unsigned int m_number_of_Chebyshev_points;

    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_position_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_rotation_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_acceleration_integrator;

    const Eigen::Matrix3d m_M_linear { m_rod_properties->getMLinear() };

    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;

};






struct DeltaInternalCouplesIntegrator : public OSNI::ODEAb {

    DeltaInternalCouplesIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                   std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                   std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;



    virtual Eigen::Vector3d computeLocalExternalCouples(unsigned int t_point) const;




    const unsigned int m_number_of_Chebyshev_points;

    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;

    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation_perturbation;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator { m_idm_integrators->m_quaternion };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_position_integrator { m_idm_integrators->m_position };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack { m_strain_parameterisation->m_K_stack };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack { m_strain_parameterisation_perturbation->m_K_stack };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack { m_strain_parameterisation->m_Lambda_stack };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack { m_strain_parameterisation_perturbation->m_Lambda_stack};



    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_rotation_integrator;


    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator { m_idm_integrators->m_angular_velocity };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_acceleration_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator { m_idm_integrators->m_linear_velocity };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;


    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator { m_idm_integrators->m_internal_forces };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_forces_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator { m_idm_integrators->m_internal_couples };


    const Eigen::Matrix3d m_M_angular { m_rod_properties->m_M.block<3,3>(0, 0) };
    const Eigen::Matrix3d m_M_linear { m_rod_properties->m_M.block<3,3>(3, 3) };




};



struct DeltaGeneralisedForcesIntegrator : public OSNI::ODEb {

    DeltaGeneralisedForcesIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_couples_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_forces_integrator;

};




struct TIDMIntegrators{

    TIDMIntegrators(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                    std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                    std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                    std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);


    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;
    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation_Delta;
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    const unsigned int m_number_of_Chebyshev_points { m_strain_parameterisation->getNumberOfChebyshewPoints() };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_rotation { std::make_shared<DeltaRotation>( m_strain_parameterisation,
                                                                                                  m_strain_parameterisation_Delta) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_position { std::make_shared<DeltaPosition>( m_strain_parameterisation,
                                                                                                  m_strain_parameterisation_Delta,
                                                                                                  m_Delta_rotation) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_angular_velocity { std::make_shared<DeltaAngularVelocity>(m_strain_parameterisation,
                                                                                                                m_strain_parameterisation_Delta,
                                                                                                                m_idm_integrators) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_linear_velocity { std::make_shared<DeltaLinearVelocity>(m_strain_parameterisation,
                                                                                                              m_strain_parameterisation_Delta,
                                                                                                              m_idm_integrators,
                                                                                                              m_Delta_angular_velocity) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_angular_acceleration { std::make_shared<DeltaAngularAccelerations>(m_strain_parameterisation,
                                                                                                                         m_strain_parameterisation_Delta,
                                                                                                                         m_idm_integrators,
                                                                                                                         m_Delta_angular_velocity) };




    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_linear_acceleration { std::make_shared<DeltaLinearAccelerations>(m_strain_parameterisation,
                                                                                                                       m_strain_parameterisation_Delta,
                                                                                                                       m_idm_integrators,
                                                                                                                       m_Delta_angular_velocity,
                                                                                                                       m_Delta_angular_acceleration,
                                                                                                                       m_Delta_linear_velocity) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_forces { std::make_shared<DeltaInternalForcesIntegrator>(m_strain_parameterisation,
                                                                                                                        m_strain_parameterisation_Delta,
                                                                                                                        m_idm_integrators,
                                                                                                                        m_rod_properties,
                                                                                                                        m_Delta_rotation,
                                                                                                                        m_Delta_angular_velocity,
                                                                                                                        m_Delta_linear_velocity,
                                                                                                                        m_Delta_linear_acceleration) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_couples { std::make_shared<DeltaInternalCouplesIntegrator>(m_strain_parameterisation,
                                                                                                                          m_strain_parameterisation_Delta,
                                                                                                                          m_idm_integrators,
                                                                                                                          m_rod_properties,
                                                                                                                          m_Delta_rotation,
                                                                                                                          m_Delta_angular_velocity,
                                                                                                                          m_Delta_angular_acceleration,
                                                                                                                          m_Delta_linear_velocity,
                                                                                                                          m_Delta_internal_forces) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_generalised_forces { std::make_shared<DeltaGeneralisedForcesIntegrator>(m_strain_parameterisation,
                                                                                                                              m_Delta_internal_couples,
                                                                                                                              m_Delta_internal_forces) };





};

}   //  namespace CROSP

#endif // TIDM_INTEGRATORS_HPP
