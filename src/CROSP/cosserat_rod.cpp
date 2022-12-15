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

#include "CROSP/CROSP/cosserat_rod.hpp"


#include <boost/math/special_functions/legendre.hpp>

namespace CROSP {

CosseratRod::CosseratRod()
{}



CosseratRod::CosseratRod(unsigned int t_number_of_Chebyshev_points)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


CosseratRod::CosseratRod(const rod_properties::RodProperties &t_rod_properties)
    : m_rod_properties( std::make_shared<rod_properties::RodProperties>(t_rod_properties) )
{}


CosseratRod::CosseratRod(const rod_properties::RodProperties &t_rod_properties,
                         unsigned int t_number_of_Chebyshev_points)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_rod_properties( std::make_shared<rod_properties::RodProperties>(t_rod_properties) )
{}


CosseratRod::CosseratRod(const polynomial_representation::PolynomialRepresentation &t_polynomial_representation)
    : m_polynomial_representation(t_polynomial_representation)
{}

CosseratRod::CosseratRod(const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                         unsigned int t_number_of_Chebyshev_points)
    : m_polynomial_representation(t_polynomial_representation),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}

CosseratRod::CosseratRod(const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                         const rod_properties::RodProperties &t_rod_properties)
    : m_polynomial_representation(t_polynomial_representation),
      m_rod_properties( std::make_shared<rod_properties::RodProperties>(t_rod_properties) )
{}

CosseratRod::CosseratRod(const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                         const rod_properties::RodProperties &t_rod_properties,
                         unsigned int t_number_of_Chebyshev_points)
    : m_polynomial_representation(t_polynomial_representation),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_rod_properties( std::make_shared<rod_properties::RodProperties>(t_rod_properties) )
{}



void CosseratRod::updateParameterisation(const Eigen::VectorXd &t_qe,
                                         const Eigen::VectorXd &t_dot_qe,
                                         const Eigen::VectorXd &t_ddot_qe)
{
    m_strain_parameterisation->updateStacks(t_qe, t_dot_qe, t_ddot_qe);
}


void CosseratRod::forwardKinematics()
{
    //  Integrate Quaternions
    m_idm_integrators->m_quaternion->solveSystem();

    //  Integrate Positions
    m_idm_integrators->m_position->solveSystem();

    //  Integrate twist
    m_idm_integrators->m_angular_velocity->solveSystem();
    m_idm_integrators->m_linear_velocity->solveSystem();

    //  Integrate accelerations
    m_idm_integrators->m_angular_acceleration->solveSystem();
    m_idm_integrators->m_linear_acceleration->solveSystem();
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

void CosseratRod::updateParameterisationVariation(const Eigen::VectorXd &t_Delta_qe,
                                                  const Eigen::VectorXd &t_Delta_dot_qe,
                                                  const Eigen::VectorXd &t_Delta_ddot_qe)
{
    m_strain_parameterisation_Delta->updateStacks(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
}

void CosseratRod::forwardTangentKinematics()
{
    //  Integrate Delta zeta
    m_tidm_integrators->m_Delta_rotation->solveSystem();
    m_tidm_integrators->m_Delta_position->solveSystem();

    //  Integrate Delta eta
    m_tidm_integrators->m_Delta_angular_velocity->solveSystem();
    m_tidm_integrators->m_Delta_linear_velocity->solveSystem();

    //  Integrate Delta dot eta
    m_tidm_integrators->m_Delta_angular_acceleration->solveSystem();
    m_tidm_integrators->m_Delta_linear_acceleration->solveSystem();
}


void CosseratRod::forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                           const Eigen::Vector3d &t_initial_Delta_position,
                                           const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                           const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                           const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                           const Eigen::Vector3d &t_initial_Delta_linear_acceleration)
{
    //  Integrate Delta zeta
    m_tidm_integrators->m_Delta_rotation->integrate( t_initial_Delta_orientation );
    m_tidm_integrators->m_Delta_position->integrate( t_initial_Delta_position );

    //  Integrate Delta eta
    m_tidm_integrators->m_Delta_angular_velocity->integrate( t_initial_Delta_angular_velocity );
    m_tidm_integrators->m_Delta_linear_velocity->integrate( t_initial_Delta_linear_velocity );

    //  Integrate Delta dot eta
    m_tidm_integrators->m_Delta_angular_acceleration->integrate( t_initial_Delta_angular_acceleration );
    m_tidm_integrators->m_Delta_linear_acceleration->integrate( t_initial_Delta_linear_acceleration );
}

::LieAlgebra::Kinematics CosseratRod::getKinematicsAtTip()const
{
    ::LieAlgebra::Kinematics rod_tip_kinematics;

    rod_tip_kinematics.m_pose = ::LieAlgebra::SE3Pose( m_idm_integrators->m_quaternion->getStateAtPoint(0),
                                                       m_idm_integrators->m_position->getStateAtPoint(0) );


    rod_tip_kinematics.m_twist << m_idm_integrators->m_angular_velocity->getStateAtPoint(0),
                                    m_idm_integrators->m_linear_velocity->getStateAtPoint(0);

    rod_tip_kinematics.m_accelerations << m_idm_integrators->m_angular_acceleration->getStateAtPoint(0),
                                            m_idm_integrators->m_linear_acceleration->getStateAtPoint(0);
    return rod_tip_kinematics;
}


void CosseratRod::backwardDynamics(const Eigen::Vector3d &t_couple_at_tip,
                                   const Eigen::Vector3d &t_force_at_tip)
{

    const auto [tip_pose, _, __] = getKinematicsAtTip();

    //  Map force and couple into local coordinates
    Eigen::Vector3d force_at_tip_local_coord = tip_pose.getRotationMatrix().transpose()*t_force_at_tip;
    Eigen::Vector3d couple_at_tip_local_coord = tip_pose.getRotationMatrix().transpose()*t_couple_at_tip;

    m_idm_integrators->m_internal_forces->integrate(force_at_tip_local_coord);
    m_idm_integrators->m_internal_couples->integrate(couple_at_tip_local_coord);


    m_idm_integrators->m_generalised_forces->integrate(Eigen::VectorXd::Zero(getCoordinatesDimension()));

}


void CosseratRod::backwardTangentDynamics(const Eigen::Vector3d &t_Delta_couple_at_tip,
                                          const Eigen::Vector3d &t_Delta_force_at_tip)
{
    const auto [tip_pose, _, __] = getKinematicsAtTip();

    //  Map force and couple into local coordinates
    Eigen::Vector3d Delta_force_at_tip_local_coord = tip_pose.getRotationMatrix().transpose()*t_Delta_force_at_tip;
    Eigen::Vector3d Delta_couple_at_tip_local_coord = tip_pose.getRotationMatrix().transpose()*t_Delta_couple_at_tip;


    m_tidm_integrators->m_Delta_internal_forces->integrate(Delta_force_at_tip_local_coord);
    m_tidm_integrators->m_Delta_internal_couples->integrate(Delta_couple_at_tip_local_coord);

    m_tidm_integrators->m_Delta_generalised_forces->integrate(Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()));
}



::LieAlgebra::Vector6d CosseratRod::IDM(const Eigen::VectorXd &t_qe,
                                        const Eigen::VectorXd &t_dot_qe,
                                        const Eigen::VectorXd &t_ddot_qe,
                                        const Eigen::Vector3d &t_couple_at_tip,
                                        const Eigen::Vector3d &t_force_at_tip)
{
    updateParameterisation(t_qe, t_dot_qe, t_ddot_qe);

    forwardKinematics();
    backwardDynamics(t_couple_at_tip,
                     t_force_at_tip);

    return getLambdaAtBase();
}

::LieAlgebra::Vector6d CosseratRod::TIDM(const Eigen::VectorXd &t_Delta_qe,
                                         const Eigen::VectorXd &t_Delta_dot_qe,
                                         const Eigen::VectorXd &t_Delta_ddot_qe,
                                         const Eigen::Vector3d &t_Delta_couple_at_tip,
                                         const Eigen::Vector3d &t_Delta_force_at_tip)
{
    updateParameterisationVariation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);


    forwardTangentKinematics();
    backwardTangentDynamics(t_Delta_couple_at_tip,
                            t_Delta_force_at_tip);

    return getDeltaLambdaAtBase();
}

Vector6d CosseratRod::getLambdaAtBase()const
{
    Vector6d Lambda;
    Lambda <<   m_idm_integrators->m_internal_couples->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN),
                m_idm_integrators->m_internal_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

    return Lambda;
}



Vector6d CosseratRod::getDeltaLambdaAtBase()const
{
    Vector6d Delta_Lambda;
    Delta_Lambda <<   m_tidm_integrators->m_Delta_internal_couples->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN),
                      m_tidm_integrators->m_Delta_internal_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

    return Delta_Lambda;
}

unsigned int CosseratRod::getCoordinatesDimension()const
{
    return m_polynomial_representation.getCoordinatesDimension();
}


Eigen::VectorXd CosseratRod::getStaticEquilibrium(const Eigen::VectorXd &t_qe) const
{

    const Eigen::VectorXd Qe = m_idm_integrators->m_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);
    return m_rod_properties->m_Kee*t_qe - Qe;
}

Eigen::VectorXd CosseratRod::getTangentStaticEquilibrium(const Eigen::VectorXd &t_Delta_qe)const
{
    return m_rod_properties->m_Kee*t_Delta_qe
            - m_tidm_integrators->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);
}





}   //  namespace CROSP
