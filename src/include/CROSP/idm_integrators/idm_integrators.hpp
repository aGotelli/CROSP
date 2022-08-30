/**
 * \file idm_integrators.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod IDM
 * \date 13-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */



#ifndef IDM_INTEGRATORS_HPP
#define IDM_INTEGRATORS_HPP

#include <memory>

#include "CROSP/rod_properties/rod_properties.hpp"
#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include "OSNI/OSNI.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/base_maps/base_maps.hpp"



namespace CROSP::idm_integrators {





struct QuaternionIntegrator : public OSNI::ODEA {

    QuaternionIntegrator(const unsigned int t_number_of_Chebyshev_points,
                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    Eigen::Matrix4d m_A_at_chebychev_point;

};



struct PositionIntegrator : public OSNI::ODEb {

    PositionIntegrator(const unsigned int t_number_of_Chebyshev_points,
                       std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    Eigen::Quaterniond m_quaternion_at_point;

};




struct AngularVelocityIntegrator : public OSNI::ODEAb {

    AngularVelocityIntegrator(const unsigned int t_number_of_Chebyshev_points,
                              std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;


};


struct LinearVelocityIntegrator : public OSNI::ODEAb {

    LinearVelocityIntegrator(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point);


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point);


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;



};




struct AngularAccelerationIntegrator : public OSNI::ODEAb {

    AngularAccelerationIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                  std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;


};







struct LinearAccelerationIntegrator : public OSNI::ODEAb {

    LinearAccelerationIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                 std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator);



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;



};




struct InternalForcesIntegrator : public OSNI::ODEAb {

    InternalForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator);



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computeDistributedForce(const unsigned int t_point) const;




    const unsigned int m_number_of_Chebyshev_points;

    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    const Eigen::Matrix3d m_M_linear { m_rod_properties->getMLinear() };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position;

};


struct InternalCouplesIntegrator : public OSNI::ODEAb {

    InternalCouplesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                              std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                              std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator);




    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computeDistributedCouple(const unsigned int t_point) const;


    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    const Eigen::Matrix3d m_M_angular { m_rod_properties->getMAngular() };

    const Eigen::Matrix3d m_M_linear { m_rod_properties->getMLinear() };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position;

    std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces;

};


struct GeneralisedForcesIntegrator : public OSNI::ODEb {

    GeneralisedForcesIntegrator(const unsigned int t_number_of_Chebyshev_points,
                                const base_maps::PolynomialRepresentation &t_polynomial_representation,
                                std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
                                std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;

};



/*!
 * \brief The IDMIntegrators struct contains all the intgrators in order to integrate the kinematics and dynamics of the rod
 *
 * This struct is a collection of integrators to use in order to compute the Forward Kinematics and the Backward Dynamics
 */
struct IDMIntegrators {

    IDMIntegrators(const unsigned int t_number_of_Chebyshev_points,
                   const base_maps::PolynomialRepresentation &t_polynomial_representation,
                   std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                   std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);

    //  Instance of the number of points to be used in the integration
    const unsigned int m_number_of_Chebyshev_points;

    //  Instance of the polynomial representation of the field of strain
    const base_maps::PolynomialRepresentation m_polynomial_representation;

    //  Instance of the strain parameterisation
    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;

    //  Instance of the rod properties
    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    //  Extract the needed matrices
    const Eigen::Matrix3d m_M_angular { m_rod_properties->getMAngular() };
    const Eigen::Matrix3d m_M_linear { m_rod_properties->getMLinear() };



    //  Integrator for the quaternions
    std::shared_ptr<OSNI::ODESolverInterface> m_quaternion { std::make_shared<QuaternionIntegrator>(m_number_of_Chebyshev_points,
                                                                                                    m_strain_parameterisation) };

    //  Integrator for the positions
    std::shared_ptr<OSNI::ODESolverInterface> m_position { std::make_shared<PositionIntegrator>(m_number_of_Chebyshev_points,
                                                                                                m_strain_parameterisation,
                                                                                                m_quaternion ) };

    //  Integrator for the angular velocities
    std::shared_ptr<OSNI::ODESolverInterface> m_angular_velocity { std::make_shared<AngularVelocityIntegrator>(m_number_of_Chebyshev_points,
                                                                                                               m_strain_parameterisation) };

    //  Integrator for the linear velocities
    std::shared_ptr<OSNI::ODESolverInterface> m_linear_velocity { std::make_shared<LinearVelocityIntegrator>(m_number_of_Chebyshev_points,
                                                                                                             m_strain_parameterisation,
                                                                                                             m_angular_velocity ) };

    //  Integrator for the angular accelerations
    std::shared_ptr<OSNI::ODESolverInterface> m_angular_acceleration { std::make_shared<AngularAccelerationIntegrator>(m_number_of_Chebyshev_points,
                                                                                                                       m_strain_parameterisation,
                                                                                                                       m_angular_velocity ) };

    //  Integrator for the linear accelerations
    std::shared_ptr<OSNI::ODESolverInterface> m_linear_acceleration { std::make_shared<LinearAccelerationIntegrator>(m_number_of_Chebyshev_points,
                                                                                                                     m_strain_parameterisation,
                                                                                                                     m_angular_velocity,
                                                                                                                     m_linear_velocity,
                                                                                                                     m_angular_acceleration ) };


    //  Integrator for the internal forces
    std::shared_ptr<OSNI::ODESolverInterface> m_internal_forces { std::make_shared<InternalForcesIntegrator>(m_number_of_Chebyshev_points,
                                                                                                             m_strain_parameterisation,
                                                                                                             m_rod_properties,
                                                                                                             m_angular_velocity,
                                                                                                             m_linear_velocity,
                                                                                                             m_linear_acceleration,
                                                                                                             m_quaternion,
                                                                                                             m_position )};

    //  Integrator for the internal couples
    std::shared_ptr<OSNI::ODESolverInterface> m_internal_couples { std::make_shared<InternalCouplesIntegrator>(m_number_of_Chebyshev_points,
                                                                                                               m_strain_parameterisation,
                                                                                                               m_rod_properties,
                                                                                                               m_angular_velocity,
                                                                                                               m_linear_velocity,
                                                                                                               m_angular_acceleration,
                                                                                                               m_quaternion,
                                                                                                               m_position,
                                                                                                               m_internal_forces )};

    //  Integrator for the generalised coordinates
    std::shared_ptr<OSNI::ODESolverInterface> m_generalised_forces { std::make_unique<GeneralisedForcesIntegrator>(m_number_of_Chebyshev_points,
                                                                                                                   m_polynomial_representation,
                                                                                                                   m_strain_parameterisation,
                                                                                                                   m_internal_couples,
                                                                                                                   m_internal_forces) };
};


}   //  namespace CROSP::IDMIntegrators





#endif // IDM_INTEGRATORS_HPP
