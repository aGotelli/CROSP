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
    m_strain_parameterisation->update(t_qe, t_dot_qe, t_ddot_qe);
}

void CosseratRod::forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                    const Eigen::Vector3d &t_initial_position,
                                    const Eigen::Vector3d &t_initial_angular_velocity,
                                    const Eigen::Vector3d &t_initial_linear_velocity,
                                    const Eigen::Vector3d &t_initial_angular_acceleration,
                                    const Eigen::Vector3d &t_initial_linear_acceleration)
{
    //  Integrate Quaternions
    m_idm_integrators->m_quaternion->integrate(t_initial_quaternion);

    //  Integrate Positions
    m_idm_integrators->m_position->integrate(t_initial_position);

    //  Integrate twist
    m_idm_integrators->m_angular_velocity->integrate(t_initial_angular_velocity);
    m_idm_integrators->m_linear_velocity->integrate(t_initial_linear_velocity);

    //  Integrate accelerations
    m_idm_integrators->m_angular_acceleration->integrate(t_initial_angular_acceleration);
    m_idm_integrators->m_linear_acceleration->integrate(t_initial_linear_acceleration);
}

void CosseratRod::updateParameterisationPerturbation(const Eigen::VectorXd &t_Delta_qe,
                                        const Eigen::VectorXd &t_Delta_dot_qe,
                                        const Eigen::VectorXd &t_Delta_ddot_qe)
{
    m_strain_parameterisation_perturbation->update(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
}

void CosseratRod::forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                           const Eigen::Vector3d &t_initial_Delta_position,
                                           const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                           const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                           const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                           const Eigen::Vector3d &t_initial_Delta_linear_acceleration)
{
    //  Integrate Delta zeta
    m_tidm_integrators->m_Delta_rotation->integrate( Eigen::Vector3d::Zero() );
    m_tidm_integrators->m_Delta_position->integrate( Eigen::Vector3d::Zero() );

    //  Integrate Delta eta
    m_tidm_integrators->m_Delta_angular_velocity->integrate(Eigen::Vector3d::Zero());
    m_tidm_integrators->m_Delta_linear_velocity->integrate(Eigen::Vector3d::Zero());

    //  Integrate Delta dot eta
    m_tidm_integrators->m_Delta_angular_acceleration->integrate(Eigen::Vector3d::Zero());
    m_tidm_integrators->m_Delta_linear_acceleration->integrate(Eigen::Vector3d::Zero());
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

}



}   //  namespace CROSP
