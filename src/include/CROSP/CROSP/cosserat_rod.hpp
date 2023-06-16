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

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"
#include "CROSP/numerical_integrators/spectral_method/spectral_integrators.hpp"
#include "CROSP/numerical_integrators/runge_kutta/runge_kutta_integrator.hpp"
#include "CROSP/numerical_integrators/magnus_expansion/magnus_integrators.hpp"

#include <boost/numeric/odeint.hpp>


namespace CROSP {


/*!
 * \brief The CosseratRod class implements the functionalities needed to simulate a Cosserat rod
 */
class CosseratRod
{
public:

    /*!
     * \brief Default constructors that initialised member with by their default initialisation
     */
    CosseratRod()=default;

    CosseratRod(polynomial_representation::PolynomialRepresentation t_polynomial_represenation,
                unsigned int t_number_of_Chebyshev_points)
        : m_polynomial_representation(t_polynomial_represenation),
          m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
    {}



    CosseratRod(polynomial_representation::PolynomialRepresentation t_polynomial_represenation,
                rod_properties::RodDimensions t_rod_dimensions,
                rod_properties::MaterialProperties t_material_properties,
                strain_parameterisation_stack::StrainFunction t_constrained_strain=strain_parameterisation_stack::default_constrained_strain);

    /*!
     * \brief updateParameterisation updates the parameterisation of the strain describing the rod shape
     * \param t_qe the set of generalised coordinates
     * \param t_dot_qe the set of the first derivatives of the generalised coordinates
     * \param t_ddot_qe the set of the second derivatives of the generalised coordinates
     *
     * This function call the StrainParameterisation::updateStacks in order to compute the current values
     * for all the K and Gamma, the angular and linear part of the strain and their derivatives
     */
    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                const Eigen::VectorXd &t_dot_qe,
                                const Eigen::VectorXd &t_ddot_qe);

    /*!
     * \brief forwardKinematics computes the forward kinematics of the rod starting from the identity pose
     *
     * This function computes the forward kinematics of the rod starting from the indentity pose and with null
     * velocities and accelerations
     */
    void forwardKinematics();


    /*!
     * \brief forwardKinematics computes the forward kinematics of the rod starting from the given states
     * \param t_initial_quaternion the initial state for the quaternion as a vector [w, x, y, z]
     * \param t_initial_position the initial state for the position
     * \param t_initial_angular_velocity the initial state for the angular velocities
     * \param t_initial_linear_velocity the initial state for the linear velocities
     * \param t_initial_angular_acceleration the initial state for the angular accelerations
     * \param t_initial_linear_acceleration the initial state for the linear accelerations
     */
    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                           const Eigen::Vector3d &t_initial_position,
                           const Eigen::Vector3d &t_initial_angular_velocity,
                           const Eigen::Vector3d &t_initial_linear_velocity,
                           const Eigen::Vector3d &t_initial_angular_acceleration,
                           const Eigen::Vector3d &t_initial_linear_acceleration);


    /*!
     * \brief forwardKinematics computes the forward kinematics of the rod starting from the given states
     * \param t_initial_quaternion the initial state for the quaternion as a Eigen::Quaterniond
     * \param t_initial_position the initial state for the position
     * \param t_initial_angular_velocity the initial state for the angular velocities
     * \param t_initial_linear_velocity the initial state for the linear velocities
     * \param t_initial_angular_acceleration the initial state for the angular accelerations
     * \param t_initial_linear_acceleration the initial state for the linear accelerations
     */
    void forwardKinematics(const Eigen::Quaterniond &t_initial_quaternion,
                           const Eigen::Vector3d &t_initial_position,
                           const Eigen::Vector3d &t_initial_angular_velocity,
                           const Eigen::Vector3d &t_initial_linear_velocity,
                           const Eigen::Vector3d &t_initial_angular_acceleration,
                           const Eigen::Vector3d &t_initial_linear_acceleration);


    /*!
     * \brief forwardKinematics computes the forward kinematics of the rod starting from the given states
     * \param t_base_pose is the pose of the first cross section in its local coordinates
     * \param t_base_twist is the twist of the first cross section in its local coordinates
     * \param t_base_acceleration is the acceleration of the first cross section in its local coordinates
     *
     * This function is a warper for forwardKinematics(const Eigen::Quaterniond &,const Eigen::Vector3d &,const Eigen::Vector3d &,const Eigen::Vector3d &,const Eigen::Vector3d &,const Eigen::Vector3d &)
     */
    inline void forwardKinematics(const ::LieAlgebra::SE3Pose &t_base_pose,
                                  const ::LieAlgebra::Vector6d &t_base_twist,
                                  const ::LieAlgebra::Vector6d &t_base_acceleration)
    {
        forwardKinematics(t_base_pose.m_quaternion, t_base_pose.m_position,
                          t_base_twist.block<3, 1>(0, 0), t_base_twist.block<3, 1>(3, 0),
                          t_base_acceleration.block<3, 1>(0, 0), t_base_acceleration.block<3, 1>(3, 0));
    }


    /*!
     * \brief updateDeltaParameterisation updates the variation of the strain due to the delta
     * \param t_Delta_qe the variation on the generalised coordinates
     * \param t_Delta_dot_qe the variation on the first derivative of the generalised coordinates
     * \param t_Delta_ddot_qe the variation on the second derivative of the  generalised coordinates
     */
    void updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                     const Eigen::VectorXd &t_Delta_dot_qe,
                                     const Eigen::VectorXd &t_Delta_ddot_qe);


    /*!
     * \brief forwardTangentKinematics computes the forward tangent kinematics of the rod starting from the null tangent state
     *
     * This function computes the forward kinematics of the rod starting from the null tangent state with null
     * delta rotation, delta position as well as velocities and accelerations
     */
    void forwardTangentKinematics();


    /*!
     * \brief forwardTangentKinematics computes the forward tangent kinematics of the rod starting from the given states
     * \param t_initial_Delta_orientation the initial variation in the orientation
     * \param t_initial_Delta_position the initial variation in the position
     * \param t_initial_Delta_angular_velocity the initial variation in the angular velocity
     * \param t_initial_Delta_linear_velocity the initial variation in the linear velocity
     * \param t_initial_Delta_angular_acceleration the initial variation in the angular acceleration
     * \param t_initial_Delta_linear_acceleration the initial variation in the linear acceleration
     */
    void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                  const Eigen::Vector3d &t_initial_Delta_position,
                                  const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                  const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                  const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                  const Eigen::Vector3d &t_initial_Delta_linear_acceleration);


    /*!
     * \brief forwardTangentKinematics computes the forward tangent kinematics of the rod starting from the given states
     * \param t_initial_Delta_pose
     * \param t_initial_Delta_twist
     * \param t_initial_Delta_acceleration
     */
    inline void forwardTangentKinematics(const ::LieAlgebra::Vector6d &t_initial_Delta_pose,
                                         const ::LieAlgebra::Vector6d &t_initial_Delta_twist,
                                         const ::LieAlgebra::Vector6d &t_initial_Delta_acceleration)
    {
        forwardTangentKinematics(t_initial_Delta_pose.block<3,1>(0, 0), t_initial_Delta_pose.block<3,1>(3, 0),
                                 t_initial_Delta_twist.block<3,1>(0, 0), t_initial_Delta_twist.block<3,1>(3, 0),
                                 t_initial_Delta_acceleration.block<3,1>(0, 0), t_initial_Delta_acceleration.block<3,1>(3, 0));
    }


    /*!
     * \brief getKinematicsAtTip gives the kinematics at the rod tip
     *
     * \return the kinematics state of the rod tip
     */
    ::LieAlgebra::Kinematics getKinematicsAtTip()const;


    /*!
     * \brief getTangentKinematicsAtTip gives the tangent kinematics at the rod tip
     *
     * \return the tangent kinematics state of the rod tip
     */
    ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const;



    /*!
     * \brief backwardDynamics computes the backward dynamics starting from the given state at the rod tip
     * \param t_Lambda_X1 is the wrench at the rod tip expressed in the local coordinates frame
     *
     * This function takes as parameter the Wrench Lambda at X=1 expressed in the frame attached to the cross section at X=1
     */
    void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1);



    /*!
     * \brief backwardTangentDynamics computes the tangent backward dynamics starting from the given state at the rod tip
     * \param t_Delta_Lambda_X1 is the tangent wrench at the rod tip expressed in the local coordinates frame
     *
     * This function takes as parameter the tangent Wrench Delta_Lambda at X=1 expressed in the frame attached to the cross section at X=1
     */
    void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1);





    /*!
     * \brief IDM computes the inverse dynamic model for the Cosserat rod
     * \param t_qe is the set of generalised coordinates
     * \param t_dot_qe is the set of the first derivative for the generalised coordinates
     * \param t_ddot_qe is the set of the second derivative for the generalised coordinates
     * \param t_initial_quaternion the initial state for the quaternion as a vector [w, x, y, z]
     * \param t_initial_position the initial position of the rod
     * \param t_initial_twist the twist at the rod base, with first the angular part followed by linear part
     * \param t_initial_acceleration the acceleration of the rod base, with first the angular part followed by linear part
     * \param t_Lambda_X1 the wrench at the rod tip in its local coordinates frame
     * \return Lambda_X0, the vector of internal wrenches expressed in the rod local frame at X=0
     */
    ::LieAlgebra::Vector6d IDM(const Eigen::VectorXd &t_qe,
                               const Eigen::VectorXd &t_dot_qe,
                               const Eigen::VectorXd &t_ddot_qe,
                               const Eigen::Vector4d &t_initial_quaternion,
                               const Eigen::Vector3d &t_initial_position,
                               const ::LieAlgebra::Vector6d &t_initial_twist,
                               const ::LieAlgebra::Vector6d &t_initial_acceleration,
                               const LieAlgebra::Vector6d &t_wrench_at_tip);



    /*!
     * \brief IDM computes the inverse dynamic model for the Cosserat rod
     * \param t_qe is the set of generalised coordinates
     * \param t_dot_qe is the set of the first derivative for the generalised coordinates
     * \param t_ddot_qe is the set of the second derivative for the generalised coordinates
     * \param t_initial_quaternion the initial state for the quaternion as a Eigen::Quaterniond
     * \param t_initial_position the initial position of the rod
     * \param t_initial_twist the twist at the rod base, with first the angular part followed by linear part
     * \param t_initial_acceleration the acceleration of the rod base, with first the angular part followed by linear part
     * \param t_Lambda_X1 the wrench at the rod tip in its local coordinates frame
     * \return Lambda_X0, the vector of internal wrenches expressed in the rod local frame at X=0
     *
     * This function is a warper around the IDM(const Eigen::VectorXd &, const Eigen::VectorXd &, const Eigen::VectorXd &, const Eigen::Vector4d &, const Eigen::Vector3d &, const ::LieAlgebra::Vector6d &, const ::LieAlgebra::Vector6d &,const LieAlgebra::Vector6d&)
     * so that it is possible to use an Eigen::Quaterniond directly in the function call.
     */
    ::LieAlgebra::Vector6d IDM(const Eigen::VectorXd &t_qe,
                               const Eigen::VectorXd &t_dot_qe,
                               const Eigen::VectorXd &t_ddot_qe,
                               const Eigen::Quaterniond &t_initial_quaternion,
                               const Eigen::Vector3d &t_initial_position,
                               const ::LieAlgebra::Vector6d &t_initial_twist,
                               const ::LieAlgebra::Vector6d &t_initial_acceleration,
                               const ::LieAlgebra::Vector6d &t_Lambda_X1);




    /*!
     * \brief TIDM computes the Tangent Inverse Dynamic Model of the rod
     * \param t_Delta_qe is the set of perturbed generalised coordinates (unitary perturbation)
     * \param t_Delta_dot_qe is the set of perdurbed derivatives of the generalised coordinates (under integrator constrains)
     * \param t_Delta_ddot_qe is the set of perdurbed double derivatives of the generalised coordinates (under integrator constrains)
     * \param t_Delta_Lambda_X1 is the unitary perturbation of the wrench at the rod tip, expressed in its own reference frame
     */
    ::LieAlgebra::Vector6d TIDM(const Eigen::VectorXd &t_Delta_qe,
                                const Eigen::VectorXd &t_Delta_dot_qe,
                                const Eigen::VectorXd &t_Delta_ddot_qe,
                                const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1);


    /*!
     * \brief getLambdaAtBase give Lambda at the rod base, expressed in local coordinates of the rod base frame
     * \return Lambda at the rod base, expressed in local coordinates of the rod base frame at X=0
     */
    ::LieAlgebra::Vector6d getLambdaAtBase()const;


    /*!
     * \brief getDeltaLambdaAtBase give Delta Lambda at the rod base, expressed in local coordinates of the rod base frame
     * \return Delta Lambda at the rod base, expressed in local coordinates of the rod base frame at X=0
     */
    ::LieAlgebra::Vector6d getDeltaLambdaAtBase()const;

    /*!
     * \brief getCoordinatesDimension gives the dimension of the rod parameterisation, namely ne*na
     * \return the dimension of the rod parameterisation, namely ne*na
     */
    inline unsigned int getCoordinatesDimension()const {return m_strain_parameterisation->m_polynomial_representation->getCoordinatesDimension();}

    /*!
     * \brief getStaticEquilibrium returns the static equilibrium of the rod Kee*qe - Q
     * \param t_qe the current set of generalised coordinates
     * \return the static equilibrium of the rod
     */
    Eigen::VectorXd getStaticInternalBalance(const Eigen::VectorXd &t_qe)const;



    Eigen::VectorXd getInternalBalance(const Eigen::VectorXd &t_qe,
                                       const Eigen::VectorXd &t_dot_qe)const;



    Eigen::VectorXd getTangentStaticInternalBalance(const Eigen::VectorXd &t_Delta_qe)const;


    Eigen::VectorXd getTangentInternalBalance(const Eigen::VectorXd &t_Delta_qe,
                                              const Eigen::VectorXd &t_Delta_dot_qe)const;


    virtual void updateInternalActuation(const double &t_current_time);

    virtual Eigen::VectorXd getQad()const;


    /*!
     * \brief getRodPositionsAtChebyshevPoints gives the position of the rod at the Chebyshev points grid
     * \return a matrix containing the stack of all the rod position at the Chebyshev points
     *
     * This function returns a stack of points. Each row contains the x, y, and z coordinates of the rod at the Chebyshev point.
     * The rows are ordered so that the first correspond to the initial position at X = 0 and the last contains the position at X = 1.
     *
     */
    inline Eigen::MatrixXd getRodPositionsAtChebyshevPoints()const{return m_cosserat_rod_integrators->getRodPositions();}


    /*!
     * \brief getRodShapeFromElasticCoordinates allows to get the shape of the rod for a given set of generalised elastic coordinates
     * \param t_qe the set of generalised elastic coordinates
     * \return a matrix containing the stack of all the rod position at the Chebyshev points
     *
     * This function makes an internal call to the function getRodPositionsAtChebyshevPoints() and shares the same return type.
     *
     */
    Eigen::MatrixXd getRodShapeFromElasticCoordinates(const Eigen::VectorXd &t_qe);


    /*!
     * \brief printProperties a function to log the rod properties with a MATLAB like layout
     */
    void printProperties();


    /*!
     * \brief updateRodLength this funtion is used to update the length of the rod
     * \param t_rod_lenght the new lenght of the rod.
     */
    [[deprecated("This function is not tested. Solve the GitHub issue before usage")]]
    virtual void updateRodLength(const double &t_rod_lenght);


    void updateRodProperties(const rod_properties::RodDimensions &t_rod_dimensions,
                             const rod_properties::MaterialProperties &t_material_properties);

    void updateRodProperties(const double &t_EI);




#ifndef DEVELOPER
protected:
#endif


    polynomial_representation::PolynomialRepresentation m_polynomial_representation;

    unsigned int m_number_of_Chebyshev_points { 21 };


    strain_parameterisation_stack::StrainFunction m_constrained_strain {
        strain_parameterisation_stack::default_constrained_strain
    };



    //  Representation of the rod via strain
    std::shared_ptr<strain_parameterisation::StrainParameterisation> m_strain_parameterisation {
        std::make_shared<strain_parameterisation::StrainParameterisation>()
    };




//#ifndef DEVELOPER
//private:
//#endif

    //  The set of rod properties
    rod_properties::RodPropertiesSPtr m_rod_properties {
        std::make_shared<rod_properties::RodProperties>(m_polynomial_representation)
    };



    //  Variables related the perturbation of the strain parameterisation
    std::shared_ptr<strain_parameterisation::StrainParameterisation> m_strain_parameterisation_Delta {
        std::make_shared<strain_parameterisation::StrainParameterisation>(m_strain_parameterisation,
                                                                          ::LieAlgebra::Vector6d::Zero())
    };


//    ::CROSP::numerical_integrators::CosseratIntegratorUPtr m_cosserat_rod_integrators {
//      std::make_unique<::CROSP::numerical_integrators::runge_kutta::RungeKuttaIntegrator>(*m_strain_parameterisation,
//                                                                                             *m_strain_parameterisation_Delta,
//                                                                                          *m_strain_parameterisation->m_polynomial_representation,
//                                                                                             m_rod_properties,
//                                                                                          m_strain_parameterisation->m_number_of_Chebyshev_points)
//    };

    ::CROSP::numerical_integrators::CosseratIntegratorUPtr m_cosserat_rod_integrators {
      std::make_unique<::CROSP::numerical_integrators::runge_kutta::RungeKuttaIntegrator>(m_polynomial_representation,
                                                                                             m_number_of_Chebyshev_points,
                                                                                             m_rod_properties)
    };



};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
