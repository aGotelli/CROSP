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


CosseratRod::CosseratRod(const unsigned int t_number_of_chebyshev_points) :
    m_strain_parameterisation( std::make_shared<strain_parameterisation::StrainParameterisation>(t_number_of_chebyshev_points) ),
    m_strain_parameterisation_perturbation( std::make_shared<strain_parameterisation::StrainParameterisationDelta>(t_number_of_chebyshev_points) )
{
    setForwardIntegratorsInitialConditions();
}


CosseratRod::CosseratRod(std::shared_ptr<material_properties::MaterialProperties> t_material_properties,
                         std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                         std::shared_ptr<strain_parameterisation::StrainParameterisationDelta> t_strain_parameterisation_perturbation) :   m_material_properties(t_material_properties),
                                                                                                            m_strain_parameterisation(t_strain_parameterisation),
                                                                                                            m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation)
{
    setForwardIntegratorsInitialConditions();
}

void CosseratRod::updateParameterisation(const Eigen::VectorXd &t_qe,
                                         const Eigen::VectorXd &t_dot_qe,
                                         const Eigen::VectorXd &t_ddot_qe)
{
    m_strain_parameterisation->update(t_qe, t_dot_qe, t_ddot_qe);
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

void CosseratRod::updateParameterisationPerturbation(const Eigen::VectorXd &t_Delta_qe,
                                        const Eigen::VectorXd &t_Delta_dot_qe,
                                        const Eigen::VectorXd &t_Delta_ddot_qe)
{
    m_strain_parameterisation_perturbation->update(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
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

    m_tidm_integrators->m_Delta_generalised_forces->integrate(Eigen::VectorXd::Zero(m_strain_parameterisation->getCoordinatesDimention()));
}


Vector6d CosseratRod::getLambdaAtBase()const
{
    Vector6d Lambda;
    Lambda <<   m_idm_integrators->m_internal_couples->getStateAtPoint(::OSNI::ROD_POSITION::BASE),
                m_idm_integrators->m_internal_forces->getStateAtPoint(::OSNI::ROD_POSITION::BASE);

    return Lambda;
}

unsigned int CosseratRod::getCoordinatesDimension()const
{
    return m_strain_parameterisation->getCoordinatesDimention();
}


Eigen::VectorXd CosseratRod::getStaticEquilibrium(const Eigen::VectorXd &t_qe) const
{
    return m_Kee*t_qe
            - m_idm_integrators->m_generalised_forces->getStateAtPoint(::OSNI::ROD_POSITION::BASE);
}


void CosseratRod::setForwardIntegratorsInitialConditions()
{
    //  Integrate Quaternions
    m_idm_integrators->m_quaternion->setInitialConditions( Eigen::Vector4d(1, 0, 0, 0) );

    //  Integrate Positions
    m_idm_integrators->m_position->setInitialConditions( Eigen::Vector3d::Zero() );

    //  Integrate twist
    m_idm_integrators->m_angular_velocity->setInitialConditions( Eigen::Vector3d::Zero() );
    m_idm_integrators->m_linear_velocity->setInitialConditions( Eigen::Vector3d::Zero() );

    //  Integrate accelerations
    m_idm_integrators->m_angular_acceleration->setInitialConditions( Eigen::Vector3d::Zero() );
    m_idm_integrators->m_linear_acceleration->setInitialConditions( Eigen::Vector3d::Zero() );

    //  Integrate Delta zeta
    m_tidm_integrators->m_Delta_rotation->setInitialConditions( Eigen::Vector3d::Zero() );
    m_tidm_integrators->m_Delta_position->setInitialConditions( Eigen::Vector3d::Zero() );

    //  Integrate Delta eta
    m_tidm_integrators->m_Delta_angular_velocity->setInitialConditions(Eigen::Vector3d::Zero());
    m_tidm_integrators->m_Delta_linear_velocity->setInitialConditions(Eigen::Vector3d::Zero());

    //  Integrate Delta dot eta
    m_tidm_integrators->m_Delta_angular_acceleration->setInitialConditions(Eigen::Vector3d::Zero());
    m_tidm_integrators->m_Delta_linear_acceleration->setInitialConditions(Eigen::Vector3d::Zero());
}

}   //  namespace CROSP
