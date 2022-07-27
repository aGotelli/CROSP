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


CosseratRod::CosseratRod(const unsigned int t_number_of_chebyshev_points) : m_number_of_chebyshev_points(t_number_of_chebyshev_points)
{}

void CosseratRod::updateParameterisation(const Eigen::VectorXd &t_qe,
                                         const Eigen::VectorXd &t_dot_qe,
                                         const Eigen::VectorXd &t_ddot_qe)
{
    Eigen::MatrixXd Phi(m_na, m_ne*m_na);

    for(unsigned int i=0; i<m_number_of_chebyshev_points; i++){
        const double point = m_Chebyshev_points[i];

        Phi = getPhi( point );

        m_K_stack->at(i) = Phi*t_qe;

        std::cout << "The point : " << i << " at : " << point << " K : " << m_K_stack->at(i).transpose() << std::endl;
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

    std::cout << "Quaternions : \n" << m_quaternion_integrator->getStack() << "\n\n\n" << std::endl;
    std::cout << "Positions : \n" << m_position_integrator->getStack() << "\n\n\n" << std::endl;
    std::cout << "Linear velocities : \n" << m_angular_velocity_integrator->getStack() << "\n\n\n" << std::endl;
    std::cout << "Angular velocities : \n" << m_linear_velocity_integrator->getStack() << "\n\n\n" << std::endl;
    std::cout << "Linear accelerations : \n" << m_angular_acceleration_integrator->getStack() << "\n\n\n" << std::endl;
    std::cout << "Angular accelerations : \n" << m_linear_acceleration_integrator->getStack() << "\n\n\n" << std::endl;

}


void CosseratRod::backwardDynamics(const Eigen::Vector3d &t_force_at_tip,
                                   const Eigen::Vector3d &t_couple_at_tip)
{
    const auto quaternion_at_tip = m_quaternion_integrator->getStateAtPoint(0);

    std::cout << "Quaternion at tip :" << quaternion_at_tip.transpose() << std::endl;


    const Eigen::Matrix3d rod_tip_orientation = Eigen::Quaterniond(m_quaternion_integrator->getStateAtPoint(0)(0),
                                                                    m_quaternion_integrator->getStateAtPoint(0)(1),
                                                                    m_quaternion_integrator->getStateAtPoint(0)(2),
                                                                    m_quaternion_integrator->getStateAtPoint(0)(3)).toRotationMatrix();
    //  Map force and couple into local coordinates
    Eigen::Vector3d force_at_tip_local_coord = rod_tip_orientation.transpose()*t_force_at_tip;

    std::cout << "Force at tip : " << t_force_at_tip.transpose() << "\n" << "Local force : " << force_at_tip_local_coord.transpose() << std::endl;

    Eigen::Vector3d couple_at_tip_local_coord = rod_tip_orientation.transpose()*t_couple_at_tip;

    std::cout << "Integrate forces \n\n";
    m_internal_forces_integrator->integrate(force_at_tip_local_coord);
    std::cout << "Integrate couples \n\n";
    m_internal_couples_integrator->integrate(couple_at_tip_local_coord);

    for(unsigned int i=1; i<5; i++)
        std::cout << "At point " << i << " force " << m_internal_forces_integrator->getStateAtPoint(i).transpose() << "\n\n";

    std::cout << "Internal forces : \n" << m_internal_forces_integrator->getStack() << "\n\n\n" << std::endl;
    std::cout << "Internal couples : \n" << m_internal_couples_integrator->getStack() << "\n\n\n" << std::endl;
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
