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

    virtual Eigen::RowVectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        m_rotation_matrix_at_point = Eigen::Quaterniond(m_quaternion_integrator->getStateAtPoint(t_point).data()).toRotationMatrix();

        return m_rotation_matrix_at_point * m_Lambda_stack->at(t_point);
    }


    std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;

    Eigen::Quaterniond m_rotation_matrix_at_point;

};

}   //  namespace CROSP





#endif // INTEGRATORS_HPP
