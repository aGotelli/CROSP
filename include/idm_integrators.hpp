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

#include "OSNI/OSNI.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"
#include "strain_parameterisation.hpp"


namespace CROSP {









struct QuaternionIntegrator : public OSNI::ODEA {

    QuaternionIntegrator(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                         const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEA(4, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD, t_number_of_Chebyshev_points),
                                                                            m_K_stack(t_K_stack)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        m_A_at_chebychev_point  <<              0             ,   -m_K_stack->at(t_point)(0),   -m_K_stack->at(t_point)(1),   -m_K_stack->at(t_point)(2),
                                    m_K_stack->at(t_point)(0) ,                0            ,    m_K_stack->at(t_point)(2),   -m_K_stack->at(t_point)(1),
                                    m_K_stack->at(t_point)(1) ,   -m_K_stack->at(t_point)(2),                0            ,    m_K_stack->at(t_point)(0),
                                    m_K_stack->at(t_point)(2) ,    m_K_stack->at(t_point)(1),   -m_K_stack->at(t_point)(0),                0            ;

        return 0.5*m_A_at_chebychev_point;
    }

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    Eigen::Matrix4d m_A_at_chebychev_point;

};



struct PositionIntegrator : public OSNI::ODEb {

    PositionIntegrator(std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                       std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                       const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD, t_number_of_Chebyshev_points),
                                                                          m_quaternion(t_quaternion_integrator),
                                                                          m_Lambda_stack(t_Lambda_stack)
    {}

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        m_quaternion_at_point = {m_quaternion->getStateAtPoint(t_point)(0),
                                m_quaternion->getStateAtPoint(t_point)(1),
                                m_quaternion->getStateAtPoint(t_point)(2),
                                m_quaternion->getStateAtPoint(t_point)(3)};

        return m_quaternion_at_point.toRotationMatrix() * m_Lambda_stack->at(t_point);
    }


    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    Eigen::Quaterniond m_quaternion_at_point;

};




struct AngularVelocityIntegrator : public OSNI::ODEAb {

    AngularVelocityIntegrator(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                              std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                              const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD, t_number_of_Chebyshev_points),
                                                                                 m_K_stack(t_K_stack),
                                                                                 m_dot_K_stack(t_dot_K_stack)
    {}


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return m_dot_K_stack->at(t_point);
    }


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;


};


struct LinearVelocityIntegrator : public OSNI::ODEAb {

    LinearVelocityIntegrator(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                             std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                             std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_Lambda_stack,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD, t_number_of_Chebyshev_points),
                                                                                m_K_stack(t_K_stack),
                                                                                m_Lambda_stack(t_Lambda_stack),
                                                                                m_dot_Lambda_stack(t_dot_Lambda_stack),
                                                                                m_angular_velocity(t_angular_velocity_integrator)
    {}


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        Eigen::Vector3d dot_Lambda = m_dot_Lambda_stack->at(t_point);
        Eigen::Matrix3d skew_Lambda = ::LieAlgebra::skew( m_Lambda_stack->at(t_point) );
        Eigen::Vector3d Omega = m_angular_velocity->getStateAtPoint(t_point);
//        return m_dot_Lambda_stack->at(t_point) - skew( m_Lambda_stack->at(t_point) )*m_angular_velocity->getStateAtPoint(t_point) ;
        return dot_Lambda - skew_Lambda*Omega;
    }


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;



};




struct AngularAccelerationIntegrator : public OSNI::ODEAb {

    AngularAccelerationIntegrator(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                                  std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                                  std::shared_ptr<const std::vector<Eigen::Vector3d>> t_ddot_K_stack,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                  const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD, t_number_of_Chebyshev_points),
                                                                                     m_K_stack(t_K_stack),
                                                                                     m_dot_K_stack(t_dot_K_stack),
                                                                                     m_ddot_K_stack(t_ddot_K_stack),
                                                                                     m_angular_velocity(t_angular_velocity_integrator)
    {}



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return m_ddot_K_stack->at(t_point) - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_angular_velocity->getStateAtPoint(t_point);
    }


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;


};







struct LinearAccelerationIntegrator : public OSNI::ODEAb {

    LinearAccelerationIntegrator(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                                 std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                                 std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                                 std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_Lambda_stack,
                                 std::shared_ptr<const std::vector<Eigen::Vector3d>> t_ddot_Lambda_stack,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                 const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD, t_number_of_Chebyshev_points),
                                                                                    m_K_stack(t_K_stack),
                                                                                    m_dot_K_stack(t_dot_K_stack),
                                                                                    m_Lambda_stack(t_Lambda_stack),
                                                                                    m_dot_Lambda_stack(t_dot_Lambda_stack),
                                                                                    m_ddot_Lambda_stack(t_ddot_Lambda_stack),
                                                                                    m_angular_velocity(t_angular_velocity_integrator),
                                                                                    m_linear_velocity(t_linear_velocity_integrator),
                                                                                    m_angular_acceleration(t_angular_acceleration_integrator)
    {}



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return m_ddot_Lambda_stack->at(t_point) - ::LieAlgebra::skew( m_Lambda_stack->at(t_point) ) * m_angular_acceleration->getStateAtPoint(t_point)
                                                - ::LieAlgebra::skew( m_dot_Lambda_stack->at(t_point) ) * m_angular_velocity->getStateAtPoint(t_point)
                                                - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_linear_velocity->getStateAtPoint(t_point);
    }


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

    InternalForcesIntegrator(const Eigen::Matrix3d t_M_linear,
                             std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                             const unsigned int t_number_of_Chebyshev_points) :                       OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD, t_number_of_Chebyshev_points),
                                                                                                      m_M_linear(t_M_linear),
                                                                                                      m_K_stack(t_K_stack),
                                                                                                      m_angular_velocity(t_angular_velocity_integrator),
                                                                                                      m_linear_velocity(t_linear_velocity_integrator),
                                                                                                      m_linear_acceleration(t_linear_acceleration_integrator),
                                                                                                      m_quaternion(t_quaternion_integrator),
                                                                                                      m_position(t_position_integrator)
    {}



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return m_M_linear*m_linear_acceleration->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_angular_velocity->getStateAtPoint(t_point) ).transpose() * m_M_linear * m_linear_velocity->getStateAtPoint(t_point)
                - computeDistributedForce(t_point);
    }


    virtual Eigen::VectorXd computeDistributedForce(const unsigned int t_point)
    {
        if(t_point == 0 or t_point == 10)
            return Eigen::Vector3d::Zero();

        Eigen::Vector4d q(m_quaternion->getStateAtPoint(t_point));

        Eigen::Matrix3d R = Eigen::Quaterniond(q(0),
                                               q(1),
                                               q(2),
                                               q(3)).toRotationMatrix();

        Eigen::Vector3d N_bar = -R.transpose()*7800*3.1416e-06*Eigen::Vector3d(0, 0, 9.81);


        return N_bar;
    }


    const Eigen::Matrix3d m_M_linear;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position;


};


struct InternalCouplesIntegrator : public OSNI::ODEAb {

    InternalCouplesIntegrator(const Eigen::Matrix3d t_M_angular,
                             const Eigen::Matrix3d t_M_linear,
                             std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                             std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                              const unsigned int t_number_of_Chebyshev_points) :                              OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD, t_number_of_Chebyshev_points),
                                                                                                              m_M_angular(t_M_angular),
                                                                                                              m_M_linear(t_M_linear),
                                                                                                              m_K_stack(t_K_stack),
                                                                                                              m_Lambda_stack(t_Lambda_stack),
                                                                                                              m_angular_velocity(t_angular_velocity_integrator),
                                                                                                              m_linear_velocity(t_linear_velocity_integrator),
                                                                                                              m_angular_acceleration(t_angular_acceleration_integrator),
                                                                                                              m_quaternion(t_quaternion_integrator),
                                                                                                              m_position(t_position_integrator),
                                                                                                              m_internal_forces(t_internal_forces_integrator)
    {}




    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return ::LieAlgebra::skew( m_Lambda_stack->at(t_point) ).transpose()*m_internal_forces->getStateAtPoint(t_point)
                + m_M_angular*m_angular_acceleration->getStateAtPoint(t_point)
                - ::LieAlgebra::skew(m_angular_velocity->getStateAtPoint(t_point)).transpose() * m_M_angular * m_angular_velocity->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_linear_velocity->getStateAtPoint(t_point) ).transpose()* m_M_linear * m_linear_velocity->getStateAtPoint(t_point)
                + computeDistributedCouple(m_quaternion->getStateAtPoint(t_point), m_position->getStateAtPoint(t_point));
    }


    virtual Eigen::VectorXd computeDistributedCouple(const Eigen::Vector4d &t_quaternion, const Eigen::Vector3d &t_position)
    {
        Eigen::Vector3d distributed_couples(0, 0, 0);
        return distributed_couples;
    }


    const Eigen::Matrix3d m_M_angular;

    const Eigen::Matrix3d m_M_linear;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position;

    std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces;


};



struct IDMIntegrators {

    IDMIntegrators(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                   const Eigen::Matrix3d &t_M_angular, const Eigen::Matrix3d &t_M_linear) : m_strain_parameterisation(t_strain_parameterisation),
                                                                                            m_M_angular(t_M_angular),
                                                                                            m_M_linear(t_M_angular)
    {}

    IDMIntegrators(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                   const Eigen::Matrix3d &t_M_angular, const Eigen::Matrix3d &t_M_linear,
                   const unsigned int t_number_of_Chebyshev_points) :   m_strain_parameterisation(t_strain_parameterisation),
                                                                        m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
                                                                        m_M_angular(t_M_angular),
                                                                        m_M_linear(t_M_angular)
    {}

    const std::shared_ptr<const StrainParameterisation> m_strain_parameterisation;
    const unsigned int m_number_of_Chebyshev_points { 17 };

    const Eigen::Matrix3d m_M_angular;
    const Eigen::Matrix3d m_M_linear;

    std::shared_ptr<OSNI::ODESolverInterface> m_quaternion { std::make_shared<QuaternionIntegrator>(m_strain_parameterisation->m_K_stack,
                                                                                                               m_number_of_Chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_position { std::make_shared<PositionIntegrator>(m_quaternion,
                                                                                                           m_strain_parameterisation->m_Lambda_stack,
                                                                                                           m_number_of_Chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_angular_velocity { std::make_shared<AngularVelocityIntegrator>(m_strain_parameterisation->m_K_stack,
                                                                                                                          m_strain_parameterisation->m_dot_K_stack,
                                                                                                                          m_number_of_Chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_linear_velocity { std::make_shared<LinearVelocityIntegrator>(m_strain_parameterisation->m_K_stack,
                                                                                                                        m_strain_parameterisation->m_Lambda_stack,
                                                                                                                        m_strain_parameterisation->m_dot_Lambda_stack,
                                                                                                                        m_angular_velocity,
                                                                                                                        m_number_of_Chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_angular_acceleration { std::make_shared<AngularAccelerationIntegrator>(m_strain_parameterisation->m_K_stack,
                                                                                                                                  m_strain_parameterisation->m_dot_K_stack,
                                                                                                                                  m_strain_parameterisation->m_ddot_K_stack,
                                                                                                                                  m_angular_velocity,
                                                                                                                                  m_number_of_Chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_linear_acceleration { std::make_shared<LinearAccelerationIntegrator>(m_strain_parameterisation->m_K_stack,
                                                                                                                                m_strain_parameterisation->m_dot_K_stack,
                                                                                                                                m_strain_parameterisation->m_Lambda_stack,
                                                                                                                                m_strain_parameterisation->m_dot_Lambda_stack,
                                                                                                                                m_strain_parameterisation->m_ddot_Lambda_stack,
                                                                                                                                m_angular_velocity,
                                                                                                                                m_linear_velocity,
                                                                                                                                m_angular_acceleration,
                                                                                                                                m_number_of_Chebyshev_points ) };


    std::shared_ptr<OSNI::ODESolverInterface> m_internal_forces { std::make_shared<InternalForcesIntegrator>(m_M_linear,
                                                                                                                        m_strain_parameterisation->m_K_stack, m_angular_velocity,
                                                                                                                        m_linear_velocity,
                                                                                                                        m_linear_acceleration,
                                                                                                                        m_quaternion,
                                                                                                                        m_position,
                                                                                                                        m_number_of_Chebyshev_points )};

    std::shared_ptr<OSNI::ODESolverInterface> m_internal_couples { std::make_shared<InternalCouplesIntegrator>(m_M_angular,
                                                                                                               m_M_linear,
                                                                                                               m_strain_parameterisation->m_K_stack,
                                                                                                               m_strain_parameterisation->m_Lambda_stack,
                                                                                                               m_angular_velocity,
                                                                                                               m_linear_velocity,
                                                                                                               m_angular_acceleration,
                                                                                                               m_quaternion,
                                                                                                               m_position,
                                                                                                               m_internal_forces,
                                                                                                               m_number_of_Chebyshev_points )};

};


}   //  namespace CROSP





#endif // IDM_INTEGRATORS_HPP
