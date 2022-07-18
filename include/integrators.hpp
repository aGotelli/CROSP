/**
 * \file integrators.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod IDM
 * \date 13-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */



#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <memory>

#include "OSNI/include/OSNI.hpp"

#include "math_tools/include/lie_algebra_utilities.hpp"


namespace CROSP {

struct QuaternionIntegrator : public OSNI::ODEA {
    QuaternionIntegrator(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack) : OSNI::ODEA(4),
                                                                                          m_K_stack(t_K_stack) {}

    QuaternionIntegrator(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                         const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEA(4, t_number_of_Chebyshev_points),
                                                                            m_K_stack(t_K_stack){}

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
                       std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack) : OSNI::ODEb(3),
                                                                                             m_quaternion_integrator(t_quaternion_integrator),
                                                                                             m_Lambda_stack(t_Lambda_stack){}

    PositionIntegrator(std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                       std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                       const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEb(3, t_number_of_Chebyshev_points),
                                                                          m_quaternion_integrator(t_quaternion_integrator),
                                                                          m_Lambda_stack(t_Lambda_stack){}

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        m_quaternion_at_point = {m_quaternion_integrator->getStateAtPoint(t_point)(0),
                                m_quaternion_integrator->getStateAtPoint(t_point)(1),
                                m_quaternion_integrator->getStateAtPoint(t_point)(2),
                                m_quaternion_integrator->getStateAtPoint(t_point)(3)};

        return m_quaternion_at_point.toRotationMatrix() * m_Lambda_stack->at(t_point);
    }


    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    Eigen::Quaterniond m_quaternion_at_point;

};




struct AngularVelocityIntegrator : public OSNI::ODEAb {

    AngularVelocityIntegrator(const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                              const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack) : OSNI::ODEAb(3),
                                                                                                         m_K_stack(t_K_stack),
                                                                                                         m_dot_K_stack(t_dot_K_stack) {}

    AngularVelocityIntegrator(const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                              const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                              const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(3, t_number_of_Chebyshev_points),
                                                                                 m_K_stack(t_K_stack),
                                                                                 m_dot_K_stack(t_dot_K_stack) {}


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -skew( m_K_stack->at(t_point) );
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return m_dot_K_stack->at(t_point);
    }


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;


};


struct LinearVelocityIntegrator : public OSNI::ODEAb {

    LinearVelocityIntegrator(const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                             const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                             const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_Lambda_stack,
                             const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator) : OSNI::ODEAb(3),
                                                                                                                    m_K_stack(t_K_stack),
                                                                                                                    m_Lambda_stack(t_Lambda_stack),
                                                                                                                    m_dot_Lambda_stack(t_dot_Lambda_stack),
                                                                                                                    m_angular_velocity_integrator(t_angular_velocity_integrator) {}

    LinearVelocityIntegrator(const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                             const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                             const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_Lambda_stack,
                             const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(3, t_number_of_Chebyshev_points),
                                                                                m_K_stack(t_K_stack),
                                                                                m_Lambda_stack(t_Lambda_stack),
                                                                                m_dot_Lambda_stack(t_dot_Lambda_stack),
                                                                                m_angular_velocity_integrator(t_angular_velocity_integrator){}


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -skew( m_K_stack->at(t_point) );
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        Eigen::Vector3d dot_Lambda = m_dot_Lambda_stack->at(t_point);
        Eigen::Matrix3d skew_Lambda = skew( m_Lambda_stack->at(t_point) );
        Eigen::Vector3d Omega = m_angular_velocity_integrator->getStateAtPoint(t_point);
//        return m_dot_Lambda_stack->at(t_point) - skew( m_Lambda_stack->at(t_point) )*m_angular_velocity_integrator->getStateAtPoint(t_point) ;
        return dot_Lambda - skew_Lambda*Omega;
    }


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;



};




struct AngularAccelerationIntegrator : public OSNI::ODEAb {

    AngularAccelerationIntegrator(const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                                  const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                                  const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_ddot_K_stack,
                                  const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator) : OSNI::ODEAb(3),
                                                                                                                         m_K_stack(t_K_stack),
                                                                                                                         m_dot_K_stack(t_dot_K_stack),
                                                                                                                         m_ddot_K_stack(t_ddot_K_stack),
                                                                                                                         m_angular_velocity_integrator(t_angular_velocity_integrator){}

    AngularAccelerationIntegrator(const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                                  const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                                  const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_ddot_K_stack,
                                  const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                  const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(3, t_number_of_Chebyshev_points),
                                                                                     m_K_stack(t_K_stack),
                                                                                     m_dot_K_stack(t_dot_K_stack),
                                                                                     m_ddot_K_stack(t_ddot_K_stack),
                                                                                     m_angular_velocity_integrator(t_angular_velocity_integrator){}



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -skew( m_K_stack->at(t_point) );
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return m_ddot_K_stack->at(t_point) - skew( m_dot_K_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point);
    }


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;


};







struct LinearAccelerationIntegrator : public OSNI::ODEAb {

    LinearAccelerationIntegrator(const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                                 const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                                 const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                                 const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_Lambda_stack,
                                 const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_ddot_Lambda_stack,
                                 const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                 const std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                 const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator) : OSNI::ODEAb(3),
                                                                                                                            m_K_stack(t_K_stack),
                                                                                                                            m_dot_K_stack(t_dot_K_stack),
                                                                                                                            m_Lambda_stack(t_Lambda_stack),
                                                                                                                            m_dot_Lambda_stack(t_dot_Lambda_stack),
                                                                                                                            m_ddot_Lambda_stack(t_ddot_Lambda_stack),
                                                                                                                            m_angular_velocity_integrator(t_angular_velocity_integrator),
                                                                                                                            m_linear_velocity_integrator(t_linear_velocity_integrator),
                                                                                                                            m_angular_acceleration_integrator(t_angular_acceleration_integrator)
    {}

    LinearAccelerationIntegrator(const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                                 const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                                 const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                                 const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_Lambda_stack,
                                 const std::shared_ptr<const std::vector<Eigen::Vector3d>> t_ddot_Lambda_stack,
                                 const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                 const std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                 const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                 const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(3, t_number_of_Chebyshev_points),
                                                                                    m_K_stack(t_K_stack),
                                                                                    m_dot_K_stack(t_dot_K_stack),
                                                                                    m_Lambda_stack(t_Lambda_stack),
                                                                                    m_dot_Lambda_stack(t_dot_Lambda_stack),
                                                                                    m_ddot_Lambda_stack(t_ddot_Lambda_stack),
                                                                                    m_angular_velocity_integrator(t_angular_velocity_integrator),
                                                                                    m_linear_velocity_integrator(t_linear_velocity_integrator),
                                                                                    m_angular_acceleration_integrator(t_angular_acceleration_integrator)
    {}



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -skew( m_K_stack->at(t_point) );
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return m_ddot_Lambda_stack->at(t_point) - skew( m_Lambda_stack->at(t_point) ) * m_angular_acceleration_integrator->getStateAtPoint(t_point)
                                                - skew( m_dot_Lambda_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point)
                                                - skew( m_dot_K_stack->at(t_point) ) * m_linear_velocity_integrator->getStateAtPoint(t_point);
    }


    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration_integrator;



};




struct InternalForcesIntegrator : public OSNI::ODEAb {
    InternalForcesIntegrator(const Eigen::Matrix3d t_M_linear,
                             std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                             const std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             const std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                             const std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator) : OSNI::ODEAb(3),
                                                                                                      m_M_linear(t_M_linear),
                                                                                                      m_K_stack(t_K_stack),
                                                                                                      m_angular_velocity_integrator(t_angular_velocity_integrator),
                                                                                                      m_linear_velocity_integrator(t_linear_velocity_integrator),
                                                                                                      m_linear_acceleration_integrator(t_linear_acceleration_integrator),
                                                                                                      m_quaternion_integrator(t_quaternion_integrator),
                                                                                                      m_position_integrator(t_position_integrator)
    {}



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {
        return -skew( m_K_stack->at(t_point) ).transpose();
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        return m_M_linear*m_linear_acceleration_integrator->getStateAtPoint(t_point)
                - skew( m_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_M_linear * m_linear_velocity_integrator->getStateAtPoint(t_point)
                + computeDistributedForce(m_quaternion_integrator->getStateAtPoint(t_point), m_position_integrator->getStateAtPoint(t_point));
    }


    virtual Eigen::VectorXd computeDistributedForce(const Eigen::Vector4d &t_quaternion, const Eigen::Vector3d &t_position)
    {
        return Eigen::Vector3d::Zero();
    }


    const Eigen::Matrix3d m_M_linear;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration_integrator;

    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator;

    std::shared_ptr<const OSNI::ODESolverInterface> m_position_integrator;


};




}   //  namespace CROSP





#endif // INTEGRATORS_HPP
