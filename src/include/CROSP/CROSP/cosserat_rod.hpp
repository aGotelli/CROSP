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



#ifndef COSSERAT_ROD_HPP
#define COSSERAT_ROD_HPP

#include <Eigen/Dense>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"

#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"

#include "CROSP/idm_integrators/idm_integrators.hpp"
#include "CROSP/tidm_integrators/tidm_integrators.hpp"


#include <boost/numeric/odeint.hpp>


namespace CROSP {


/*!
 * \brief The CosseratRod class implements the functionalities needed to simulate a Cosserat rod
 */
class CosseratRod
{
public:

    CosseratRod()=default;


    CosseratRod(const strain_parameterisation::StrainParameterisation &t_strain_parameterisation,
                const rod_properties::RodProperties &t_rod_properties);







//    /*!
//     * \brief updateParameterisation updates the parameterisation of the strain describing the rod shape
//     * \param t_qe the set of generalised coordinates
//     * \param t_dot_qe the set of the first derivatives of the generalised coordinates
//     * \param t_ddot_qe the set of the second derivatives of the generalised coordinates
//     *
//     * This function call the StrainParameterisation::updateStacks in order to compute the current values
//     * for all the K and Gamma, the angular and linear part of the strain and their derivatives
//     */
//    void updateParameterisation(const Eigen::VectorXd &t_qe,
//                                const Eigen::VectorXd &t_dot_qe,
//                                const Eigen::VectorXd &t_ddot_qe);

//    /*!
//     * \brief forwardKinematics computes the forward kinematics of the rod starting from the identity pose
//     *
//     * This function computes the forward kinematics of the rod starting from the indentity pose and with null
//     * velocities and accelerations
//     */
//    void forwardKinematics();


//    /*!
//     * \brief forwardKinematics computes the forward kinematics of the rod starting from the given states
//     * \param t_initial_quaternion the initial state for the quaternion as a vector [w, x, y, z]
//     * \param t_initial_position the initial state for the position
//     * \param t_initial_angular_velocity the initial state for the angular velocities
//     * \param t_initial_linear_velocity the initial state for the linear velocities
//     * \param t_initial_angular_acceleration the initial state for the angular accelerations
//     * \param t_initial_linear_acceleration the initial state for the linear accelerations
//     */
//    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
//                           const Eigen::Vector3d &t_initial_position,
//                           const Eigen::Vector3d &t_initial_angular_velocity,
//                           const Eigen::Vector3d &t_initial_linear_velocity,
//                           const Eigen::Vector3d &t_initial_angular_acceleration,
//                           const Eigen::Vector3d &t_initial_linear_acceleration);


//    /*!
//     * \brief forwardKinematics computes the forward kinematics of the rod starting from the given states
//     * \param t_initial_quaternion the initial state for the quaternion as a Eigen::Quaterniond
//     * \param t_initial_position the initial state for the position
//     * \param t_initial_angular_velocity the initial state for the angular velocities
//     * \param t_initial_linear_velocity the initial state for the linear velocities
//     * \param t_initial_angular_acceleration the initial state for the angular accelerations
//     * \param t_initial_linear_acceleration the initial state for the linear accelerations
//     */
//    void forwardKinematics(const Eigen::Quaterniond &t_initial_quaternion,
//                           const Eigen::Vector3d &t_initial_position,
//                           const Eigen::Vector3d &t_initial_angular_velocity,
//                           const Eigen::Vector3d &t_initial_linear_velocity,
//                           const Eigen::Vector3d &t_initial_angular_acceleration,
//                           const Eigen::Vector3d &t_initial_linear_acceleration);


//    /*!
//     * \brief forwardKinematics computes the forward kinematics of the rod starting from the given states
//     * \param t_base_pose is the pose of the first cross section in its local coordinates
//     * \param t_base_twist is the twist of the first cross section in its local coordinates
//     * \param t_base_acceleration is the acceleration of the first cross section in its local coordinates
//     *
//     * This function is a warper for forwardKinematics(const Eigen::Quaterniond &,const Eigen::Vector3d &,const Eigen::Vector3d &,const Eigen::Vector3d &,const Eigen::Vector3d &,const Eigen::Vector3d &)
//     */
//    inline void forwardKinematics(const ::LieAlgebra::SE3Pose &t_base_pose,
//                                  const ::LieAlgebra::Vector6d &t_base_twist,
//                                  const ::LieAlgebra::Vector6d &t_base_acceleration)
//    {
//        forwardKinematics(t_base_pose.m_quaternion, t_base_pose.m_position,
//                          t_base_twist.block<3, 1>(0, 0), t_base_twist.block<3, 1>(3, 0),
//                          t_base_acceleration.block<3, 1>(0, 0), t_base_acceleration.block<3, 1>(3, 0));
//    }


//    /*!
//     * \brief updateParameterisationVariation updates the variation of the strain due to the delta
//     * \param t_Delta_qe the variation on the generalised coordinates
//     * \param t_Delta_dot_qe the variation on the first derivative of the generalised coordinates
//     * \param t_Delta_ddot_qe the variation on the second derivative of the  generalised coordinates
//     */
//    void updateParameterisationVariation(const Eigen::VectorXd &t_Delta_qe,
//                                         const Eigen::VectorXd &t_Delta_dot_qe,
//                                         const Eigen::VectorXd &t_Delta_ddot_qe);


//    /*!
//     * \brief forwardTangentKinematics computes the forward tangent kinematics of the rod starting from the null tangent state
//     *
//     * This function computes the forward kinematics of the rod starting from the null tangent state with null
//     * delta rotation, delta position as well as velocities and accelerations
//     */
//    void forwardTangentKinematics();


//    /*!
//     * \brief forwardTangentKinematics computes the forward tangent kinematics of the rod starting from the given states
//     * \param t_initial_Delta_orientation the initial variation in the orientation
//     * \param t_initial_Delta_position the initial variation in the position
//     * \param t_initial_Delta_angular_velocity the initial variation in the angular velocity
//     * \param t_initial_Delta_linear_velocity the initial variation in the linear velocity
//     * \param t_initial_Delta_angular_acceleration the initial variation in the angular acceleration
//     * \param t_initial_Delta_linear_acceleration the initial variation in the linear acceleration
//     */
//    void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
//                                  const Eigen::Vector3d &t_initial_Delta_position,
//                                  const Eigen::Vector3d &t_initial_Delta_angular_velocity,
//                                  const Eigen::Vector3d &t_initial_Delta_linear_velocity,
//                                  const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
//                                  const Eigen::Vector3d &t_initial_Delta_linear_acceleration);


//    /*!
//     * \brief forwardTangentKinematics computes the forward tangent kinematics of the rod starting from the given states
//     * \param t_initial_Delta_pose
//     * \param t_initial_Delta_twist
//     * \param t_initial_Delta_acceleration
//     */
//    inline void forwardTangentKinematics(const ::LieAlgebra::Vector6d &t_initial_Delta_pose,
//                                         const ::LieAlgebra::Vector6d &t_initial_Delta_twist,
//                                         const ::LieAlgebra::Vector6d &t_initial_Delta_acceleration)
//    {
//        forwardTangentKinematics(t_initial_Delta_pose.block<3,1>(0, 0), t_initial_Delta_pose.block<3,1>(3, 0),
//                                 t_initial_Delta_twist.block<3,1>(0, 0), t_initial_Delta_twist.block<3,1>(3, 0),
//                                 t_initial_Delta_acceleration.block<3,1>(0, 0), t_initial_Delta_acceleration.block<3,1>(3, 0));
//    }


//    /*!
//     * \brief getKinematicsAtTip gives the kinematics at the rod tip
//     *
//     * \return the kinematics state of the rod tip
//     */
//    ::LieAlgebra::Kinematics getKinematicsAtTip()const;


//    /*!
//     * \brief getTangentKinematicsAtTip gives the tangent kinematics at the rod tip
//     *
//     * \return the tangent kinematics state of the rod tip
//     */
//    ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const;


//    /*!
//     * \brief backwardDynamics computes the backward dynamics starting from the given state at the rod tip
//     * \param t_couple_at_tip The couple expressed in reference coordinates that acts on the rod tip
//     * \param t_force_at_tip The force expressed in reference coordinates that acts on the rod tip
//     */
//    [[deprecated("This function is deprecated as its usage is not clear.\n It is not clear for the user how to deal with projections and the frame of reference")]]
//    void backwardDynamics(const Eigen::Vector3d &t_couple_at_tip,
//                          const Eigen::Vector3d &t_force_at_tip);

//    /*!
//     * \brief backwardDynamics computes the backward dynamics starting from the given state at the rod tip
//     * \param t_Lambda_X1 is the wrench at the rod tip expressed in the local coordinates frame
//     *
//     * This function takes as parameter the Wrench Lambda at X=1 expressed in the frame attached to the cross section at X=1
//     */
//    void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1);


//    /*!
//     * \brief backwardTangentDynamics computes the tangent backward dynamics starting from the given state at the rod tip
//     * \param t_Delta_couple_at_tip The tangent couple expressed in reference coordinates that acts on the rod tip
//     * \param t_Delta_force_at_tip The tangent force expressed in reference coordinates that acts on the rod tip
//     */
//    [[deprecated("This function is deprecated as its usage is not clear.\n It is not clear for the user how to deal with projections and the frame of reference")]]
//    void backwardTangentDynamics(const Eigen::Vector3d &t_Delta_couple_at_tip,
//                                 const Eigen::Vector3d &t_Delta_force_at_tip);

//    /*!
//     * \brief backwardTangentDynamics computes the tangent backward dynamics starting from the given state at the rod tip
//     * \param t_Delta_Lambda_X1 is the tangent wrench at the rod tip expressed in the local coordinates frame
//     *
//     * This function takes as parameter the tangent Wrench Delta_Lambda at X=1 expressed in the frame attached to the cross section at X=1
//     */
//    void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1);


//    ::LieAlgebra::Vector6d IDM(const Eigen::VectorXd &t_qe,
//                               const Eigen::VectorXd &t_dot_qe,
//                               const Eigen::VectorXd &t_ddot_qe,
//                               const Eigen::Vector3d &t_couple_at_tip,
//                               const Eigen::Vector3d &t_force_at_tip);


//    /*!
//     * \brief IDM computes the inverse dynamic model for the Cosserat rod
//     * \param t_qe is the set of generalised coordinates
//     * \param t_dot_qe is the set of the first derivative for the generalised coordinates
//     * \param t_ddot_qe is the set of the second derivative for the generalised coordinates
//     * \param t_initial_quaternion the initial state for the quaternion as a vector [w, x, y, z]
//     * \param t_initial_position the initial position of the rod
//     * \param t_initial_twist the twist at the rod base, with first the angular part followed by linear part
//     * \param t_initial_acceleration the acceleration of the rod base, with first the angular part followed by linear part
//     * \param t_wrench_at_tip the wrench at the rod tip, with the couple first and then the forces
//     * \return Lambda, the vector of internal wrenches expressed in the rod local frame
//     */
//    ::LieAlgebra::Vector6d IDM(const Eigen::VectorXd &t_qe,
//                               const Eigen::VectorXd &t_dot_qe,
//                               const Eigen::VectorXd &t_ddot_qe,
//                               const Eigen::Vector4d &t_initial_quaternion,
//                               const Eigen::Vector3d &t_initial_position,
//                               const ::LieAlgebra::Vector6d &t_initial_twist,
//                               const ::LieAlgebra::Vector6d &t_initial_acceleration,
//                               const LieAlgebra::Vector6d &t_wrench_at_tip);



//    /*!
//     * \brief IDM computes the inverse dynamic model for the Cosserat rod
//     * \param t_qe is the set of generalised coordinates
//     * \param t_dot_qe is the set of the first derivative for the generalised coordinates
//     * \param t_ddot_qe is the set of the second derivative for the generalised coordinates
//     * \param t_initial_quaternion the initial state for the quaternion as a Eigen::Quaterniond
//     * \param t_initial_position the initial position of the rod
//     * \param t_initial_twist the twist at the rod base, with first the angular part followed by linear part
//     * \param t_initial_acceleration the acceleration of the rod base, with first the angular part followed by linear part
//     * \param t_wrench_at_tip the wrench at the rod tip, with the couple first and then the forces
//     * \return Lambda, the vector of internal wrenches expressed in the rod local frame
//     */
//    ::LieAlgebra::Vector6d IDM(const Eigen::VectorXd &t_qe,
//                               const Eigen::VectorXd &t_dot_qe,
//                               const Eigen::VectorXd &t_ddot_qe,
//                               const Eigen::Quaterniond &t_initial_quaternion,
//                               const Eigen::Vector3d &t_initial_position,
//                               const ::LieAlgebra::Vector6d &t_initial_twist,
//                               const ::LieAlgebra::Vector6d &t_initial_acceleration,
//                               const LieAlgebra::Vector6d &t_wrench_at_tip);



//    ::LieAlgebra::Vector6d TIDM(const Eigen::VectorXd &t_Delta_qe,
//                                const Eigen::VectorXd &t_Delta_dot_qe,
//                                const Eigen::VectorXd &t_Delta_ddot_qe,
//                                const Eigen::Vector3d &t_Delta_couple_at_tip,
//                                const Eigen::Vector3d &t_Delta_force_at_tip);


//    /*!
//     * \brief getLambdaAtBase give Lambda at the rod base, expressed in local coordinates of the rod base frame
//     * \return Lambda at the rod base, expressed in local coordinates of the rod base frame
//     */
//    ::LieAlgebra::Vector6d getLambdaAtBase()const;


//    /*!
//     * \brief getDeltaLambdaAtBase give Delta Lambda at the rod base, expressed in local coordinates of the rod base frame
//     * \return Delta Lambda at the rod base, expressed in local coordinates of the rod base frame
//     */
//    ::LieAlgebra::Vector6d getDeltaLambdaAtBase()const;

//    /*!
//     * \brief getCoordinatesDimension gives the dimension of the rod parameterisation, namely ne*na
//     * \return the dimension of the rod parameterisation, namely ne*na
//     */
//    unsigned int getCoordinatesDimension()const;

//    /*!
//     * \brief getStaticEquilibrium returns the static equilibrium of the rod Kee*qe - Q
//     * \param t_qe the current set of generalised coordinates
//     * \return the static equilibrium of the rod
//     */
//    Eigen::VectorXd getStaticEquilibrium(const Eigen::VectorXd &t_qe)const;


//    Eigen::VectorXd getStaticEquilibrium(const Eigen::VectorXd &t_qe,
//                                         const Eigen::VectorXd &t_dot_qe) const
//    {
//        Eigen::MatrixXd Dee = m_rod_properties->m_Dee;
//        Eigen::MatrixXd Kee = m_rod_properties->m_Kee;


//        Eigen::VectorXd Qe = Kee * t_qe;
//        Eigen::VectorXd Ce = Dee * t_dot_qe;
//        Eigen::VectorXd Qa = m_idm_integrators->m_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);
//        return Qe + Ce - Qa;

////        return m_rod_properties->m_Kee*t_qe
////                + m_rod_properties->m_Dee*t_dot_qe
////                - m_idm_integrators->m_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);
//    }



//    Eigen::VectorXd getTangentStaticEquilibrium(const Eigen::VectorXd &t_Delta_qe,
//                                                const Eigen::VectorXd &t_Delta_dot_qe)const
//    {
//        Eigen::MatrixXd Dee = m_rod_properties->m_Dee;
//        Eigen::MatrixXd Kee = m_rod_properties->m_Kee;


//        Eigen::VectorXd Delta_Qe = Kee * t_Delta_qe;
//        Eigen::VectorXd Delta_Ce = Dee * t_Delta_dot_qe;
//        Eigen::VectorXd Delta_Qa = m_tidm_integrators->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);
//        return Delta_Qe + Delta_Ce - Delta_Qa;

////        return m_rod_properties->m_Kee*t_Delta_qe
////                + m_rod_properties->m_Dee*t_Delta_dot_qe
////                - m_tidm_integrators->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);
//    }


//    Eigen::VectorXd getTangentStaticEquilibrium(const Eigen::VectorXd &t_Delta_qe)const;



//    inline Eigen::MatrixXd getRodPositionsAtChebyshevPoints()const{return m_idm_integrators->m_position->getStackAsMatrix();}


//    Eigen::MatrixXd getRodShape(const Eigen::VectorXd &t_qe)const;


//    void printProperties();

#ifndef DEVELOPER
private:
#endif


    //  Representation of the rod via strain
    const std::shared_ptr<strain_parameterisation::StrainParameterisation> m_strain_parameterisation;



    //  The set of rod properties
    std::shared_ptr<rod_properties::RodProperties> m_rod_properties;

//    //  All the strain releted variables
////    std::shared_ptr<strain_parameterisation::StrainParameterisation> m_strain_parameterisation { [this](){

////        auto strain_parameterisation = std::make_shared<strain_parameterisation::StrainParameterisation>(m_polynomial_representation,
////                                                                                                         m_number_of_Chebyshev_points);

////        const unsigned int ne = m_polynomial_representation.getCoordinatesDimension();
////        strain_parameterisation->updateStacks(Eigen::VectorXd::Zero(ne),
////                                              Eigen::VectorXd::Zero(ne),
////                                              Eigen::VectorXd::Zero(ne));

////        return strain_parameterisation;
////    }() };


//    //  Variables related the perturbation of the strain parameterisation
//    std::shared_ptr<strain_parameterisation::StrainParameterisation> m_strain_parameterisation_Delta {
//        std::make_shared<strain_parameterisation::StrainParameterisation>(m_polynomial_representation,
//                                                                          ::LieAlgebra::Vector6d::Zero(),
//                                                                          m_number_of_Chebyshev_points) };
//    //  The set of integrators needed for the IDM
//    std::shared_ptr<idm_integrators::IDMIntegrators> m_idm_integrators {
//        std::make_shared<idm_integrators::IDMIntegrators>(m_number_of_Chebyshev_points,
//                                                          m_polynomial_representation,
//                                                          m_strain_parameterisation,
//                                                          m_rod_properties )};

//    //  The set of integrators needed for the TIDM
//    std::shared_ptr<tidm_integrators::TIDMIntegrators> m_tidm_integrators {
//        std::make_shared<tidm_integrators::TIDMIntegrators>(m_number_of_Chebyshev_points,
//                                                            m_polynomial_representation,
//                                                            m_strain_parameterisation,
//                                                            m_strain_parameterisation_Delta,
//                                                            m_idm_integrators,
//                                                            m_rod_properties) };



};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
