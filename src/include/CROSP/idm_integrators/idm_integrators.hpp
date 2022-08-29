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




namespace CROSP::idm_integrators {





struct QuaternionIntegrator : public OSNI::ODEA {

    QuaternionIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation);

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    Eigen::Matrix4d m_A_at_chebychev_point;

};



struct PositionIntegrator : public OSNI::ODEb {

    PositionIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    Eigen::Quaterniond m_quaternion_at_point;

};




struct AngularVelocityIntegrator : public OSNI::ODEAb {

    AngularVelocityIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;


};


struct LinearVelocityIntegrator : public OSNI::ODEAb {

    LinearVelocityIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point);


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point);


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;



};




struct AngularAccelerationIntegrator : public OSNI::ODEAb {

    AngularAccelerationIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator);


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;


};







struct LinearAccelerationIntegrator : public OSNI::ODEAb {

    LinearAccelerationIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
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

    InternalForcesIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
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

    InternalCouplesIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
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

    GeneralisedForcesIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
                                std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                                const unsigned int t_number_of_Chebyshev_points);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;

};



struct IDMIntegrators {

    IDMIntegrators(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                   std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);



    const std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;


    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;
    const Eigen::Matrix3d m_M_angular { m_rod_properties->getMAngular() };
    const Eigen::Matrix3d m_M_linear { m_rod_properties->getMLinear() };

    const unsigned int m_number_of_Chebyshev_points { m_strain_parameterisation->getNumberOfChebyshewPoints() };

    std::shared_ptr<OSNI::ODESolverInterface> m_quaternion { std::make_shared<QuaternionIntegrator>( m_strain_parameterisation ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_position { std::make_shared<PositionIntegrator>(m_strain_parameterisation,
                                                                                                m_quaternion ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_angular_velocity { std::make_shared<AngularVelocityIntegrator>(m_strain_parameterisation) };

    std::shared_ptr<OSNI::ODESolverInterface> m_linear_velocity { std::make_shared<LinearVelocityIntegrator>(m_strain_parameterisation,
                                                                                                             m_angular_velocity ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_angular_acceleration { std::make_shared<AngularAccelerationIntegrator>(m_strain_parameterisation,
                                                                                                                       m_angular_velocity ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_linear_acceleration { std::make_shared<LinearAccelerationIntegrator>(m_strain_parameterisation,
                                                                                                                     m_angular_velocity,
                                                                                                                     m_linear_velocity,
                                                                                                                     m_angular_acceleration ) };


    std::shared_ptr<OSNI::ODESolverInterface> m_internal_forces { std::make_shared<InternalForcesIntegrator>(m_strain_parameterisation,
                                                                                                             m_rod_properties,
                                                                                                             m_angular_velocity,
                                                                                                             m_linear_velocity,
                                                                                                             m_linear_acceleration,
                                                                                                             m_quaternion,
                                                                                                             m_position )};

    std::shared_ptr<OSNI::ODESolverInterface> m_internal_couples { std::make_shared<InternalCouplesIntegrator>(m_strain_parameterisation,
                                                                                                               m_rod_properties,
                                                                                                               m_angular_velocity,
                                                                                                               m_linear_velocity,
                                                                                                               m_angular_acceleration,
                                                                                                               m_quaternion,
                                                                                                               m_position,
                                                                                                               m_internal_forces )};


    std::shared_ptr<OSNI::ODESolverInterface> m_generalised_forces { std::make_unique<GeneralisedForcesIntegrator>(m_strain_parameterisation,
                                                                                                                   m_internal_couples,
                                                                                                                   m_internal_forces,
                                                                                                                   m_number_of_Chebyshev_points) };
};


}   //  namespace CROSP::IDMIntegrators





#endif // IDM_INTEGRATORS_HPP
