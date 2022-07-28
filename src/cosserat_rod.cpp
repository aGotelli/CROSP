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


CosseratRod::CosseratRod(const unsigned int t_number_of_chebyshev_points) : m_strain_parameterisation( std::make_shared<StrainParameterisation>(t_number_of_chebyshev_points) )
{}

void CosseratRod::updateParameterisation(const Eigen::VectorXd &t_qe,
                                         const Eigen::VectorXd &t_dot_qe,
                                         const Eigen::VectorXd &t_ddot_qe)
{
//    Eigen::MatrixXd Phi(m_na, m_ne*m_na);

//    for(unsigned int i=0; i<m_number_of_chebyshev_points; i++){
//        const double point = m_Chebyshev_points[i];

//        Phi = getPhi( point );

//        m_K_stack->at(i) = Phi*t_qe;

//        std::cout << "The point : " << i << " at : " << point << " K : " << m_K_stack->at(i).transpose() << std::endl;
//        m_dot_K_stack->at(i) = Phi*t_dot_qe;
//        m_ddot_K_stack->at(i) = Phi*t_ddot_qe;

//    }

//    for(unsigned int i=0; i<m_number_of_chebyshev_points; i++){
//        m_Lambda_stack->at(i) = Eigen::Vector3d(1, 0, 0);
//        m_dot_Lambda_stack->at(i) = Eigen::Vector3d::Zero();
//        m_ddot_Lambda_stack->at(i) = Eigen::Vector3d::Zero();
//    }


    m_strain_parameterisation->update(t_qe, t_dot_qe, t_ddot_qe);
}

void CosseratRod::forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                    const Eigen::Vector3d &t_initial_position,
                                    const Eigen::Vector3d &t_initial_angular_velocity,
                                    const Eigen::Vector3d &t_initial_linear_velocity,
                                    const Eigen::Vector3d &t_initial_angular_acceleration,
                                    const Eigen::Vector3d &t_initial_linear_acceleration)
{
    m_idm_integrators->m_quaternion->integrate(t_initial_quaternion);

    m_idm_integrators->m_position->integrate(t_initial_position);

    m_idm_integrators->m_angular_velocity->integrate(t_initial_angular_velocity);

    m_idm_integrators->m_linear_velocity->integrate(t_initial_linear_velocity);

    m_idm_integrators->m_angular_acceleration->integrate(t_initial_angular_acceleration);

    m_idm_integrators->m_linear_acceleration->integrate(t_initial_linear_acceleration);

//    std::cout << "Quaternions : \n" << m_idm_integrators->m_quaternion->getStack() << "\n\n\n" << std::endl;
//    std::cout << "Positions : \n" << m_idm_integrators->m_position->getStack() << "\n\n\n" << std::endl;
//    std::cout << "Linear velocities : \n" << m_idm_integrators->m_angular_velocity->getStack() << "\n\n\n" << std::endl;
//    std::cout << "Angular velocities : \n" << m_idm_integrators->m_linear_velocity->getStack() << "\n\n\n" << std::endl;
//    std::cout << "Linear accelerations : \n" << m_idm_integrators->m_angular_acceleration->getStack() << "\n\n\n" << std::endl;
//    std::cout << "Angular accelerations : \n" << m_idm_integrators->m_linear_acceleration->getStack() << "\n\n\n" << std::endl;

}


::LieAlgebra::Kinematics CosseratRod::getKinematicsAtTip()const
{
    ::LieAlgebra::Kinematics rod_tip_kinematics;

    rod_tip_kinematics.m_pose = ::LieAlgebra::SE3Pose( m_idm_integrators->m_quaternion->getStateAtPoint(0),
                                                       m_idm_integrators->m_position->getStateAtPoint(0) );


    rod_tip_kinematics.m_twist = ::LieAlgebra::Screw( m_idm_integrators->m_angular_velocity->getStateAtPoint(0),
                                                      m_idm_integrators->m_linear_velocity->getStateAtPoint(0) );

    rod_tip_kinematics.m_accelerations = ::LieAlgebra::Screw( m_idm_integrators->m_angular_acceleration->getStateAtPoint(0),
                                                              m_idm_integrators->m_linear_acceleration->getStateAtPoint(0) );
    return rod_tip_kinematics;
}


void CosseratRod::backwardDynamics(const Eigen::Vector3d &t_force_at_tip,
                                   const Eigen::Vector3d &t_couple_at_tip)
{

    const auto [tip_pose, _, __] = getKinematicsAtTip();

    //  Map force and couple into local coordinates
    Eigen::Vector3d force_at_tip_local_coord = tip_pose.getRotationMatrix().transpose()*t_force_at_tip;


    Eigen::Vector3d couple_at_tip_local_coord = tip_pose.getRotationMatrix().transpose()*t_couple_at_tip;

    m_idm_integrators->m_internal_forces->integrate(force_at_tip_local_coord);
    m_idm_integrators->m_internal_couples->integrate(couple_at_tip_local_coord);

//    std::cout << "Internal forces : \n" << m_idm_integrators->m_internal_forces->getStack() << "\n\n\n" << std::endl;
//    std::cout << "Internal couples : \n" << m_idm_integrators->m_internal_couples->getStack() << "\n\n\n" << std::endl;
}

//Eigen::MatrixXd CosseratRod::getPhi(const double& t_X,
//                                    const double& t_begin,
//                                    const double& t_end)const
//{
//    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
//    double x = ( 2 * t_X - ( t_end + t_begin) ) / ( t_end - t_begin );

//    //  Compute the values of the polynomial for every element of the strain field
//    Eigen::VectorXd Phi_i(m_ne, 1);
//    for(unsigned int i=0; i<m_ne; i++)
//        Phi_i[i] = m_polynomial_base(i, x);


//    //  Define the matrix of bases
//    Eigen::MatrixXd Phi = Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_na, m_na), Phi_i.transpose());


//    return Phi;
//}


}   //  namespace CROSP
