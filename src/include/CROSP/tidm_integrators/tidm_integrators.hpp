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
    DeltaRotation(const unsigned int t_number_of_Chebyshev_points,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta);

    DeltaRotation(const unsigned int t_number_of_Chebyshev_points,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                  const Eigen::Vector3d t_initial_condition);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
};


struct DeltaPosition : public OSNI::ODEAb {

    DeltaPosition(const unsigned int t_number_of_Chebyshev_points,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                  std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation);

    DeltaPosition(const unsigned int t_number_of_Chebyshev_points,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                  std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation,
                  const Eigen::Vector3d t_initial_condition);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_delta_rotation;
};




struct DeltaAngularVelocity : public OSNI::ODEAb {
    DeltaAngularVelocity(const unsigned int t_number_of_Chebyshev_points,
                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                         std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators);

    DeltaAngularVelocity(const unsigned int t_number_of_Chebyshev_points,
                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                         std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                         const Eigen::Vector3d t_initial_condition);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;

};


struct DeltaLinearVelocity : public OSNI::ODEAb {
    DeltaLinearVelocity(const unsigned int t_number_of_Chebyshev_points,
                        std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                        std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                        std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                        std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator);

    DeltaLinearVelocity(const unsigned int t_number_of_Chebyshev_points,
                        std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                        std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                        std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                        std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                        const Eigen::Vector3d t_initial_condition);

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

    DeltaAngularAccelerations(const unsigned int t_number_of_Chebyshev_points,
                              std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                              std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                              std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator);

    DeltaAngularAccelerations(const unsigned int t_number_of_Chebyshev_points,
                              std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                              std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                              std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                              const Eigen::Vector3d t_initial_condition);

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

    DeltaLinearAccelerations(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                             std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator);

    DeltaLinearAccelerations(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                             std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                             const Eigen::Vector3d t_initial_condition);

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

    DeltaInternalForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                  std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                  std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator);

    DeltaInternalForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                  std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                  std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator,
                                  const Eigen::Vector3d t_initial_condition);

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

    DeltaInternalCouplesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                   std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                   std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator);

    DeltaInternalCouplesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_perturbation,
                                   std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                                   std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                   const Eigen::Vector3d t_initial_condition);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;



    virtual Eigen::Vector3d computeLocalExternalCouples(unsigned int) const;




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

    DeltaGeneralisedForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                     const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator);

    DeltaGeneralisedForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                     const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                                     std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                     const Eigen::Vector3d t_initial_condition);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    Eigen::MatrixXd m_B;

    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_couples_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_forces_integrator;

};




struct TIDMIntegrators{

    TIDMIntegrators(const unsigned int t_number_of_Chebyshev_points,
                    const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                    std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                    std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                    std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                    std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);

    TIDMIntegrators(const unsigned int t_index,
                    const unsigned int t_number_of_Chebyshev_points,
                    const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                    std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                    std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                    std::shared_ptr<const idm_integrators::IDMIntegrators> t_idm_integrators,
                    std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);


    const unsigned int m_index {0};

    //  Instance of the number of points to be used in the integration
    const unsigned int m_number_of_Chebyshev_points;

    //  Instance of the polynomial representation of the field of strain
    const polynomial_representation::PolynomialRepresentation m_polynomial_representation;

    //  Instance of the strain parameterisation
    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;

    //  Instance of the variation of the strain parameterisation
    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation_Delta;

    //  Instance of the IDM integrators
    const std::shared_ptr<const idm_integrators::IDMIntegrators> m_idm_integrators;

    //  Instance of the rod properties
    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;



    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_rotation { std::make_shared<DeltaRotation>(m_number_of_Chebyshev_points,
                                                                                                 m_strain_parameterisation,
                                                                                                 m_strain_parameterisation_Delta,
                                                                                                 Eigen::Vector3d::Zero()) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_position { std::make_shared<DeltaPosition>(m_number_of_Chebyshev_points,
                                                                                                 m_strain_parameterisation,
                                                                                                 m_strain_parameterisation_Delta,
                                                                                                 m_Delta_rotation,
                                                                                                 Eigen::Vector3d::Zero()) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_angular_velocity { std::make_shared<DeltaAngularVelocity>(m_number_of_Chebyshev_points,
                                                                                                                m_strain_parameterisation,
                                                                                                                m_strain_parameterisation_Delta,
                                                                                                                m_idm_integrators,
                                                                                                                Eigen::Vector3d::Zero()) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_linear_velocity { std::make_shared<DeltaLinearVelocity>(m_number_of_Chebyshev_points,
                                                                                                              m_strain_parameterisation,
                                                                                                              m_strain_parameterisation_Delta,
                                                                                                              m_idm_integrators,
                                                                                                              m_Delta_angular_velocity,
                                                                                                              Eigen::Vector3d::Zero()) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_angular_acceleration { std::make_shared<DeltaAngularAccelerations>(m_number_of_Chebyshev_points,
                                                                                                                         m_strain_parameterisation,
                                                                                                                         m_strain_parameterisation_Delta,
                                                                                                                         m_idm_integrators,
                                                                                                                         m_Delta_angular_velocity,
                                                                                                                         Eigen::Vector3d::Zero()) };




    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_linear_acceleration { std::make_shared<DeltaLinearAccelerations>(m_number_of_Chebyshev_points,
                                                                                                                       m_strain_parameterisation,
                                                                                                                       m_strain_parameterisation_Delta,
                                                                                                                       m_idm_integrators,
                                                                                                                       m_Delta_angular_velocity,
                                                                                                                       m_Delta_angular_acceleration,
                                                                                                                       m_Delta_linear_velocity,
                                                                                                                       Eigen::Vector3d::Zero()) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_forces { std::make_shared<DeltaInternalForcesIntegrator>(m_number_of_Chebyshev_points,
                                                                                                                        m_strain_parameterisation,
                                                                                                                        m_strain_parameterisation_Delta,
                                                                                                                        m_idm_integrators,
                                                                                                                        m_rod_properties,
                                                                                                                        m_Delta_rotation,
                                                                                                                        m_Delta_angular_velocity,
                                                                                                                        m_Delta_linear_velocity,
                                                                                                                        m_Delta_linear_acceleration,
                                                                                                                        Eigen::Vector3d::Zero()) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_couples { std::make_shared<DeltaInternalCouplesIntegrator>(m_number_of_Chebyshev_points,
                                                                                                                          m_strain_parameterisation,
                                                                                                                          m_strain_parameterisation_Delta,
                                                                                                                          m_idm_integrators,
                                                                                                                          m_rod_properties,
                                                                                                                          m_Delta_rotation,
                                                                                                                          m_Delta_angular_velocity,
                                                                                                                          m_Delta_angular_acceleration,
                                                                                                                          m_Delta_linear_velocity,
                                                                                                                          m_Delta_internal_forces,
                                                                                                                          Eigen::Vector3d::Zero()) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_generalised_forces { std::make_shared<DeltaGeneralisedForcesIntegrator>(m_number_of_Chebyshev_points,
                                                                                                                              m_polynomial_representation,
                                                                                                                              m_strain_parameterisation,
                                                                                                                              m_Delta_internal_couples,
                                                                                                                              m_Delta_internal_forces,
                                                                                                                              Eigen::Vector3d::Zero()) };





};

}   //  namespace CROSP

#endif // TIDM_INTEGRATORS_HPP
