/**
 * \file cosserat_rod.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the definition of strain based parameterised Cosserat rod
 * \date 12-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */

#include "cosserat_rod.hpp"


#include <boost/math/special_functions/legendre.hpp>

namespace CROSP {

CosseratRod::CosseratRod() :
    m_Chebyshev_points( ComputeChebyshevPoints(m_number_of_chebyshev_points) ),
    m_polynomial_base([](const unsigned int t_point, const double& t_x) {
                            return boost::math::legendre_p(t_point, t_x);
                        }),
    m_K_stack( std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) ),
    m_Lambda_stack( std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) ),
    m_dot_K_stack( std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) ),
    m_dot_Lambda_stack( std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) ),
    m_ddot_K_stack( std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) ),
    m_ddot_Lambda_stack( std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) ),
    //  The integrators
    m_quaternion_integrator(std::make_shared<QuaternionIntegrator>(m_K_stack)),
    m_position_integrator(std::make_shared<PositionIntegrator>(m_quaternion_integrator, m_Lambda_stack)),
    m_angular_velocity_integrator(std::make_shared<AngularVelocityIntegrator>(m_K_stack, m_dot_K_stack)),
    m_linear_velocity_integrator(std::make_shared<LinearVelocityIntegrator>(m_K_stack, m_Lambda_stack, m_dot_Lambda_stack, m_angular_velocity_integrator)),
    m_angular_acceleration_integrator(std::make_shared<AngularAccelerationIntegrator>(m_K_stack, m_dot_K_stack, m_ddot_K_stack, m_angular_velocity_integrator)),
    m_linear_acceleration_integrator(std::make_shared<LinearAccelerationIntegrator>(m_K_stack, m_dot_K_stack, m_Lambda_stack, m_dot_Lambda_stack, m_ddot_Lambda_stack, m_angular_velocity_integrator,
                                                                                     m_linear_velocity_integrator, m_angular_acceleration_integrator))
{}

void CosseratRod::updateParameterisation(const Eigen::VectorXd &t_qe,
                                         const Eigen::VectorXd &t_dot_qe,
                                         const Eigen::VectorXd &t_ddot_qe)
{
    Eigen::MatrixXd Phi(m_na, m_ne*m_na);

    for(unsigned int i=0; i<m_number_of_chebyshev_points; i++){
        Phi = getPhi(m_Chebyshev_points[i]);

        m_K_stack->at(i) = Phi*t_qe;
        m_dot_K_stack->at(i) = Phi*t_dot_qe;
        m_ddot_K_stack->at(i) = Phi*t_ddot_qe;

    }

    for(unsigned int i=0; i<m_number_of_chebyshev_points; i++){
        m_Lambda_stack->at(i) = Eigen::Vector3d(1, 0, 0);
        m_dot_Lambda_stack->at(i) = Eigen::Vector3d::Zero();
        m_ddot_Lambda_stack->at(i) = Eigen::Vector3d::Zero();
    }
}

void CosseratRod::forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                    const Eigen::Vector3d &t_initial_position,
                                    const Eigen::Vector3d &t_initial_angular_velocity,
                                    const Eigen::Vector3d &t_initial_linear_velocity,
                                    const Eigen::Vector3d &t_initial_angular_acceleration,
                                    const Eigen::Vector3d &t_initial_linear_acceleration)
{
    m_quaternion_integrator->integrate(t_initial_quaternion);

    m_position_integrator->integrate(t_initial_position);

    m_angular_velocity_integrator->integrate(t_initial_angular_velocity);

    m_linear_velocity_integrator->integrate(t_initial_linear_velocity);

    m_angular_acceleration_integrator->integrate(t_initial_angular_acceleration);

    m_linear_acceleration_integrator->integrate(t_initial_linear_acceleration);

}


Eigen::MatrixXd CosseratRod::getPhi(const double& t_X,
                                    const double& t_begin,
                                    const double& t_end) const
{
    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    double x = ( 2 * t_X - ( t_end + t_begin) ) / ( t_end - t_begin );

    //  Compute the values of the polynomial for every element of the strain field
    Eigen::VectorXd Phi_i(m_ne, 1);
    for(unsigned int i=0; i<m_ne; i++)
        Phi_i[i] = m_polynomial_base(i, x);


    //  Define the matrix of bases
    Eigen::MatrixXd Phi = Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_na, m_na), Phi_i.transpose());


    return Phi;
}



}   //  namespace CROSP
