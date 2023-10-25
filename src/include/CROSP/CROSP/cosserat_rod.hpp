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

//#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"
#include "CROSP/numerical_integrators/spectral_method/spectral_integrators.hpp"
#include "CROSP/numerical_integrators/explicit_methods/explicit_methods.hpp"
#include "CROSP/numerical_integrators/magnus_expansion/magnus_integrators.hpp"

#include <boost/numeric/odeint.hpp>


namespace CROSP {


typedef ::CROSP::numerical_integrators::spectral_method::SpectralIntegrators Spectral;
typedef ::CROSP::numerical_integrators::explicit_methods::ExplicitIntegrator<> Ode45;
typedef ::CROSP::numerical_integrators::magnus_expansion::MagnusIntegrators Magnus;


/*!
 * \brief The CosseratRod class implements the functionalities needed to simulate a Cosserat rod
 *
 * By default, this class uses the spectral integration of SpectralIntegrators
 */
template<numerical_integrators::CosseratIntegrator NumericalIntegrator=Spectral>
class CosseratRod
{
public:

    /*!
     * \brief Default constructors that initialised member with by their default initialisation
     */
    CosseratRod()=default;




    CosseratRod(polynomial_representation::PolynomialRepresentation t_polynomial_represenation,
                unsigned int t_number_of_Chebyshev_points,
                rod_properties::RodDimensions t_rod_dimensions=rod_properties::RodDimensions(),
                rod_properties::MaterialProperties t_material_properties=rod_properties::MaterialProperties(),
                strain_parameterisation_stack::StrainFunction t_constrained_strain=strain_parameterisation_stack::default_constrained_strain)
        : m_polynomial_representation(t_polynomial_represenation),
          m_rod_properties(
              std::make_shared<rod_properties::RodProperties>(t_rod_dimensions,
                                                              t_material_properties)
              ),
          m_cosserat_rod_integrators (
                std::make_unique<NumericalIntegrator>(m_polynomial_representation,
                                                     t_number_of_Chebyshev_points,
                                                     m_rod_properties,
                                                     t_constrained_strain)
              )
    {}


    /*!
     * \brief addInternalActuation adds an internal actuation to the Cosserat rod
     * \param t_distributed_actuation the internal distributed actuation acting on the rod
     *
     * This function adds an internal distributed actuation to the rod.
     * It modifies the integrators routines accordingly in order to obtain the projection of the
     * internal actuation on the generalised coordinates space; namely Q_ad.
     * In case this function is not called, Q_ad=0 and does not add any step to the numerical integrations
     */
    void addInternalActuation(::ATORS::distributed_actuation::DistributedActuationUptr t_distributed_actuation)
    {
        m_cosserat_rod_integrators->addInternalActuation( std::move(t_distributed_actuation) );
    }


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
                                const Eigen::VectorXd &t_ddot_qe)
    {
        m_cosserat_rod_integrators->updateParameterisation(t_qe, t_dot_qe, t_ddot_qe);
    }

    /*!
     * \brief forwardKinematics computes the forward kinematics of the rod starting from the identity pose
     *
     * This function computes the forward kinematics of the rod starting from the indentity pose and with null
     * velocities and accelerations
     */
    void forwardKinematics()
    {
        m_cosserat_rod_integrators->forwardKinematics();
    }


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
                           const Eigen::Vector3d &t_initial_linear_acceleration)
    {
        m_cosserat_rod_integrators->forwardKinematics(t_initial_quaternion,
                                                      t_initial_position,
                                                      t_initial_angular_velocity,
                                                      t_initial_linear_velocity,
                                                      t_initial_angular_acceleration,
                                                      t_initial_linear_acceleration);
    }



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
                           const Eigen::Vector3d &t_initial_linear_acceleration)
    {
        const Eigen::Vector4d Q(t_initial_quaternion.w(),
                                t_initial_quaternion.x(),
                                t_initial_quaternion.y(),
                                t_initial_quaternion.z());

        forwardKinematics(Q, t_initial_position,
                          t_initial_angular_velocity, t_initial_linear_velocity,
                          t_initial_angular_acceleration, t_initial_linear_acceleration);
    }


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
                                     const Eigen::VectorXd &t_Delta_ddot_qe)
    {
        m_cosserat_rod_integrators->updateDeltaParameterisation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
    }


    /*!
     * \brief forwardTangentKinematics computes the forward tangent kinematics of the rod starting from the null tangent state
     *
     * This function computes the forward kinematics of the rod starting from the null tangent state with null
     * delta rotation, delta position as well as velocities and accelerations
     */
    void forwardTangentKinematics()
    {
        m_cosserat_rod_integrators->forwardTangentKinematics();
    }


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
                                  const Eigen::Vector3d &t_initial_Delta_linear_acceleration)
    {
        m_cosserat_rod_integrators->forwardTangentKinematics(t_initial_Delta_orientation,
                                                             t_initial_Delta_position,
                                                             t_initial_Delta_angular_velocity,
                                                             t_initial_Delta_linear_velocity,
                                                             t_initial_Delta_angular_acceleration,
                                                             t_initial_Delta_linear_acceleration);
    }


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
    ::LieAlgebra::Kinematics getKinematicsAtTip()const
    {
        return m_cosserat_rod_integrators->getKinematicsAtTip();
    }


    /*!
     * \brief getTangentKinematicsAtTip gives the tangent kinematics at the rod tip
     *
     * \return the tangent kinematics state of the rod tip
     */
    ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const
    {
        return m_cosserat_rod_integrators->getTangentKinematicsAtTip();
    }



    /*!
     * \brief backwardDynamics computes the backward dynamics starting from the given state at the rod tip
     * \param t_Lambda_X1 is the wrench at the rod tip expressed in the local coordinates frame
     *
     * This function takes as parameter the Wrench Lambda at X=1 expressed in the frame attached to the cross section at X=1
     */
    void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1)
    {
        m_cosserat_rod_integrators->backwardDynamics(t_Lambda_X1);
    }



    /*!
     * \brief backwardTangentDynamics computes the tangent backward dynamics starting from the given state at the rod tip
     * \param t_Delta_Lambda_X1 is the tangent wrench at the rod tip expressed in the local coordinates frame
     *
     * This function takes as parameter the tangent Wrench Delta_Lambda at X=1 expressed in the frame attached to the cross section at X=1
     */
    void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
    {
        m_cosserat_rod_integrators->backwardTangentDynamics(t_Delta_Lambda_X1);

    }





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
                               const LieAlgebra::Vector6d &t_Lambda_X1)
    {
        m_cosserat_rod_integrators->updateParameterisation(t_qe, t_dot_qe, t_ddot_qe);

        m_cosserat_rod_integrators->forwardKinematics(t_initial_quaternion,
                                                      t_initial_position,
                                                      t_initial_twist.block<3, 1>(0, 0),
                                                      t_initial_twist.block<3, 1>(3, 0),
                                                      t_initial_acceleration.block<3, 1>(0, 0),
                                                      t_initial_acceleration.block<3, 1>(3, 0));

        m_cosserat_rod_integrators->backwardDynamics(t_Lambda_X1);

        return getLambdaAtBase();
    }




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
                               const ::LieAlgebra::Vector6d &t_Lambda_X1)
    {
        const Eigen::Vector4d initial_quaternion(t_initial_quaternion.w(),
                                                 t_initial_quaternion.x(),
                                                 t_initial_quaternion.y(),
                                                 t_initial_quaternion.z());

        return IDM(t_qe,
                   t_dot_qe,
                   t_ddot_qe,
                   initial_quaternion,
                   t_initial_position,
                   t_initial_twist,
                   t_initial_acceleration,
                   t_Lambda_X1);
    }




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
                                const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
    {
        updateDeltaParameterisation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);


        forwardTangentKinematics();
        backwardTangentDynamics(t_Delta_Lambda_X1);

        return getDeltaLambdaAtBase();
    }


    /*!
     * \brief getLambdaAtBase give Lambda at the rod base, expressed in local coordinates of the rod base frame
     * \return Lambda at the rod base, expressed in local coordinates of the rod base frame at X=0
     */
    ::LieAlgebra::Vector6d getLambdaAtBase()const
    {
        return m_cosserat_rod_integrators->getLambdaAtBase();
    }


    /*!
     * \brief getDeltaLambdaAtBase give Delta Lambda at the rod base, expressed in local coordinates of the rod base frame
     * \return Delta Lambda at the rod base, expressed in local coordinates of the rod base frame at X=0
     */
    ::LieAlgebra::Vector6d getDeltaLambdaAtBase()const
    {
        return m_cosserat_rod_integrators->getDeltaLambdaAtBase();
    }

    /*!
     * \brief getCoordinatesDimension gives the dimension of the rod parameterisation, namely ne*na
     * \return the dimension of the rod parameterisation, namely ne*na
     */
    inline unsigned int getCoordinatesDimension()const
    {
        return m_polynomial_representation.getCoordinatesDimension();
    }

    /*!
     * \brief getStaticEquilibrium returns the static equilibrium of the rod Kee*qe - Q
     * \param t_qe the current set of generalised coordinates
     * \return the static equilibrium of the rod
     */
    Eigen::VectorXd getStaticInternalBalance(const Eigen::VectorXd &t_qe)const
    {

        const Eigen::VectorXd Qe = m_Kee * t_qe;
        const Eigen::VectorXd Qa = m_cosserat_rod_integrators->getQaAtBase();
        const Eigen::VectorXd Q_ad = m_cosserat_rod_integrators->getQad();

        const Eigen::VectorXd internal_balance = Qe - Qa - Q_ad;
        return internal_balance;
    }



    Eigen::VectorXd getInternalBalance(const Eigen::VectorXd &t_qe,
                                       const Eigen::VectorXd &t_dot_qe)const
    {

        const Eigen::VectorXd Qe = m_Dee * t_qe;
        const Eigen::VectorXd Ce = m_Kee * t_dot_qe;
        const Eigen::VectorXd Qa = m_cosserat_rod_integrators->getQaAtBase();
        const Eigen::VectorXd Q_ad = m_cosserat_rod_integrators->getQad();

        const Eigen::VectorXd internal_balance = Qe + Ce - Qa - Q_ad;
        return internal_balance;

    }



    Eigen::VectorXd getTangentStaticInternalBalance(const Eigen::VectorXd &t_Delta_qe)const
    {

        Eigen::VectorXd Delta_Qe = m_Kee * t_Delta_qe;
        Eigen::VectorXd Delta_Qa = m_cosserat_rod_integrators->getDeltaQaAtBase();

        Eigen::VectorXd Delta_internal_balance = Delta_Qe - Delta_Qa;

        return Delta_internal_balance;
    }


    Eigen::VectorXd getTangentInternalBalance(const Eigen::VectorXd &t_Delta_qe,
                                              const Eigen::VectorXd &t_Delta_dot_qe)const
    {

        Eigen::VectorXd Delta_Qe = m_Kee * t_Delta_qe;
        Eigen::VectorXd Delta_Ce = m_Dee * t_Delta_dot_qe;
        Eigen::VectorXd Delta_Qa = m_cosserat_rod_integrators->getDeltaQaAtBase();


        Eigen::VectorXd Delta_internal_balance = Delta_Qe + Delta_Ce - Delta_Qa;
        return Delta_internal_balance;
    }


    void updateInternalActuation(const double &t_current_time)
    {
        m_cosserat_rod_integrators->updateInternalActuation(t_current_time);
    }



    /*!
     * \brief getRodPositionsAtChebyshevPoints gives the position of the rod at the Chebyshev points grid
     * \return a matrix containing the stack of all the rod position at the Chebyshev points
     *
     * This function returns a stack of points. Each row contains the x, y, and z coordinates of the rod at the Chebyshev point.
     * The rows are ordered so that the first correspond to the initial position at X = 0 and the last contains the position at X = 1.
     *
     */
    inline Eigen::MatrixXd getRodPositionsAtChebyshevPoints()const
    {
        return m_cosserat_rod_integrators->getRodPositions();
    }


    /*!
     * \brief getRodShapeFromElasticCoordinates allows to get the shape of the rod for a given set of generalised elastic coordinates
     * \param t_qe the set of generalised elastic coordinates
     * \return a matrix containing the stack of all the rod position at the Chebyshev points
     *
     * This function makes an internal call to the function getRodPositionsAtChebyshevPoints() and shares the same return type.
     *
     */
    Eigen::MatrixXd getRodShapeFromElasticCoordinates(const Eigen::VectorXd &t_qe)
    {

        updateParameterisation(t_qe, 0*t_qe, 0*t_qe);

        forwardKinematics();

        return getRodPositionsAtChebyshevPoints();
    }


    /*!
     * \brief printProperties a function to log the rod properties with a MATLAB like layout
     */
    void printProperties()
    {
        std::stringstream rod_properties;



        {   //  Start listing rod dimensions
        const auto l = m_rod_properties->m_rod_dimensions.m_L;
        rod_properties << "Rod dimensions :\n"
                          "     l : " << l << "\n" <<
                          m_rod_properties->m_rod_dimensions.m_cross_section->printProperties();
        }

        {   //  Start listing material properties
        const auto E = m_rod_properties->m_material_properties.m_E;
        const auto G = m_rod_properties->m_material_properties.m_G;
        const auto rho = m_rod_properties->m_material_properties.m_rho;
        rod_properties << "Material properties :\n"
                          "     E : " << std::setprecision(2) << std::scientific << E << "\n"
                          "     G : " << std::setprecision(2) << std::scientific << G << "\n"
                          "     ρ : " << std::setprecision(0) << std::fixed << rho << "\n";
        }



        {   //  Strain parameterisation
        rod_properties << "Strain parameterisation : \n";

        Eigen::MatrixXi parameterization_details_matrix = Eigen::MatrixXi::Zero(2, 6);


        Eigen::VectorXi def(6);
        for(unsigned int i=0; i<6; i++)
            parameterization_details_matrix.row(0)[i] = m_polynomial_representation.m_admitted_deformations[i];


        unsigned int j =0;
        for(unsigned int i=0; i<m_polynomial_representation.m_admitted_deformations.size(); i++)
            if(m_polynomial_representation.m_admitted_deformations[i])
                parameterization_details_matrix.row(1)[i] = m_polynomial_representation.m_number_of_modes_stack[j++];


        std::stringstream parameterization_details;
        parameterization_details << parameterization_details_matrix;

        std::string log = parameterization_details.str();
        auto newline = log.find('\n');
        std::string admitted_deformations = log.substr(0, newline);
        std::string number_of_modes = log.substr(++newline);


            //  Print the admitted deformations
        rod_properties << "     Rod deformations : " << admitted_deformations << "\n";


           //  Also print the number of modes
        rod_properties << "     Number of modes  : " << number_of_modes << "\n";



           //  Details about the base
        const auto poly_base = m_polynomial_representation.m_polynomial_base;

        rod_properties << "     Polynomial base : ";

        std::string base = "Custom base";
        const unsigned int p = 5;
        const double x = 0.5;
        if( poly_base(p, x) == ::CROSP::polynomial_representation::legendre_polynomial_base(p, x) )
            base = "Legendre";

        if( poly_base(p, x) == ::CROSP::polynomial_representation::chebyshev_polynomial_base(p, x) )
            base = "Chebyshev";

        rod_properties << base << "\n";


        }



        {   //  Integrators


            rod_properties << "Numerical Integration:\n" << m_cosserat_rod_integrators->printIntegratorProperties();

        }



        std::cout << rod_properties.str() << "\n";
        std::cout.flush();
    }


    /*!
     * \brief updateRodLength this funtion is used to update the length of the rod
     * \param t_rod_lenght the new lenght of the rod.
     */
    [[deprecated("This function is not tested. Solve the GitHub issue before usage")]]
    void updateRodLength(const double &t_rod_lenght)
    {
        const double scale = t_rod_lenght/m_rod_properties->m_rod_dimensions.m_L;
        m_rod_properties->m_rod_dimensions.m_L = t_rod_lenght;

        m_Kee *= scale;
        m_Dee *= scale;

        m_cosserat_rod_integrators->updateIntegrationDomain(t_rod_lenght);
    }


    [[deprecated("This function is not tested. Solve the GitHub issue before usage")]]
    void updateRodProperties(const rod_properties::RodDimensions &t_rod_dimensions,
                             const rod_properties::MaterialProperties &t_material_properties)
    {
        Eigen::MatrixXd Ha = m_polynomial_representation.m_B.transpose() * m_rod_properties->m_H * m_polynomial_representation.m_B;
        double L = m_rod_properties->m_rod_dimensions.m_L;

        const auto scale_matrix = L*Ha;


        m_rod_properties =
            std::make_shared<rod_properties::RodProperties>(t_rod_dimensions,
                                                            t_material_properties);


        Ha = m_polynomial_representation.m_B.transpose() * m_rod_properties->m_H * m_polynomial_representation.m_B;
        L = m_rod_properties->m_rod_dimensions.m_L;

        const auto ratio_matrix = (L*Ha).inverse()*scale_matrix;


        unsigned int row_index = 0;
        unsigned int col_index = 0;
        for(unsigned int i=0; i<m_polynomial_representation.m_na; i++){

            unsigned int ne_i = m_polynomial_representation.m_number_of_modes_stack[i];

            m_Kee.block(row_index, col_index, ne_i, ne_i) *= ratio_matrix(i, i);
        }

        m_Dee = m_rod_properties->m_material_properties.m_mu*m_Kee;

        m_cosserat_rod_integrators->updateIntegrationDomain( L );


    }

//    void updateRodProperties(const double &t_EI);




#ifndef DEVELOPER
protected:
#endif


    polynomial_representation::PolynomialRepresentation m_polynomial_representation;


    //  The set of rod properties
    rod_properties::RodPropertiesSPtr m_rod_properties;



    std::unique_ptr<NumericalIntegrator> m_cosserat_rod_integrators;



    /// \brief m_Kee The generalised elasticity matrix
    Eigen::MatrixXd m_Kee { [this]()->Eigen::MatrixXd
        {

            const auto B = m_polynomial_representation.m_B;
            const Eigen::MatrixXd Ha = B.transpose() * m_rod_properties->m_H * B;


            Eigen::MatrixXd Kee = m_cosserat_rod_integrators->integratePhiTPhi();



            const auto L = m_rod_properties->m_rod_dimensions.m_L;

            unsigned int Kee_index = 0;
            for(unsigned int i=0; i<m_polynomial_representation.m_na; i++){
                double scale = L*Ha(i, i);

                unsigned int ne_i = m_polynomial_representation.m_number_of_modes_stack[i];

                Kee.block(Kee_index, Kee_index, ne_i, ne_i) *= scale;

                Kee_index += ne_i;
            }


            return Kee;
        }()
    };

    /// \brief m_Dee The matrix of the internal dumping
    Eigen::MatrixXd m_Dee { m_rod_properties->m_material_properties.m_mu*m_Kee };



};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
