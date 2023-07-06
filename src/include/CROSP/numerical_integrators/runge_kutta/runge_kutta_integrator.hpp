#ifndef RUNGEKUTTAINTEGRATOR_HPP
#define RUNGEKUTTAINTEGRATOR_HPP



#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"


namespace CROSP::numerical_integrators::runge_kutta {


/*!
 * \brief The RungeKuttaIntegrator class provides a numeric integrator for the Cosserat IDM and TIDM using the standard Runge-Kutta integrator
 */
class RungeKuttaIntegrator
{
public:
    RungeKuttaIntegrator(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                         const unsigned int t_number_of_Chebyshev_points,
                         std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);


    std::string printIntegratorProperties()const;


    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                        const Eigen::VectorXd &t_dot_qe,
                                        const Eigen::VectorXd &t_ddot_qe);


    void forwardKinematics();


    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                   const Eigen::Vector3d &t_initial_position,
                                   const Eigen::Vector3d &t_initial_angular_velocity,
                                   const Eigen::Vector3d &t_initial_linear_velocity,
                                   const Eigen::Vector3d &t_initial_angular_acceleration,
                                   const Eigen::Vector3d &t_initial_linear_acceleration);


    void updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                             const Eigen::VectorXd &t_Delta_dot_qe,
                                             const Eigen::VectorXd &t_Delta_ddot_qe);


    void forwardTangentKinematics();


    void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                          const Eigen::Vector3d &t_initial_Delta_position,
                                          const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                          const Eigen::Vector3d &t_initial_Delta_linear_acceleration);


    ::LieAlgebra::Kinematics getKinematicsAtTip()const;


    ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const;


    void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1);


    void updateInternalActuation(const double &t_time);


    Eigen::VectorXd getQad() const;


    void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1);


    ::LieAlgebra::Vector6d getLambdaAtBase()const;


    LieAlgebra::Vector6d getDeltaLambdaAtBase()const;


    ::LieAlgebra::Vector6d getQaAtBase()const;


    LieAlgebra::Vector6d getDeltaQaAtBase()const;


    void updateIntegrationDomain(const double &t_rod_lenght);


    Eigen::MatrixXd getRodPositions()const;


    FullODEStatesObservations getIDMStatesObservations();

    FullODEStatesObservations getTIDMStatesObservations();


private:


    ::LieAlgebra::Vector6d m_constant_strain { ::LieAlgebra::Vector6d::Unit(3) };

    polynomial_representation::PolynomialRepresentation m_polynomial_representation;
    std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    double m_dX { 0.005 };

    unsigned int m_ne { m_polynomial_representation.getCoordinatesDimension() };


    double m_rod_length { m_rod_properties->m_rod_dimensions.m_L };


    unsigned int m_number_of_Chebyshev_points { 31 };
    std::vector<double> m_Chebyshev_points_reversed { [this](){
            std::vector<double> Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points);

            std::vector<double> Chebyshev_points_reversed;


            for(int i=m_number_of_Chebyshev_points-1; i>=0; i--)
                Chebyshev_points_reversed.push_back( Chebyshev_points[i] );

            return Chebyshev_points_reversed;
        }()
    };








    Eigen::VectorXd m_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };
    Eigen::VectorXd m_dot_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };
    Eigen::VectorXd m_ddot_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };


    Eigen::VectorXd m_Delta_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };
    Eigen::VectorXd m_Delta_dot_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };
    Eigen::VectorXd m_Delta_ddot_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };




    //  Define a compiled time known dimension for the state of the pose (quaternion + position)
    typedef Eigen::Matrix<double, 7, 1> PoseState;

    typedef Eigen::Matrix<double, 19, 1> ForwardKinematicState;

    ForwardKinematicState m_forward_integration_state_X0 { ForwardKinematicState::Zero() };
    ForwardKinematicState m_forward_integration_state_X1 { ForwardKinematicState::Zero() };



    Eigen::VectorXd m_backward_integration_state_X0 { Eigen::VectorXd::Zero(25+m_ne) };
    Eigen::VectorXd m_backward_integration_state_X1 { Eigen::VectorXd::Zero(25+m_ne) };


    PoseState forwardStaticODEs(const RungeKuttaIntegrator::PoseState &t_state,
                                const ::LieAlgebra::Vector6d &t_Xi) const;



    ForwardKinematicState forwardODEs(const ForwardKinematicState &t_state,
                                      const ::LieAlgebra::Vector6d &t_Xi,
                                      const ::LieAlgebra::Vector6d &t_dot_Xi,
                                      const ::LieAlgebra::Vector6d &t_ddot_Xi) const;



    ::LieAlgebra::Vector6d getLambdaPrime(const Eigen::Quaterniond &t_Q,
                                          const ::LieAlgebra::Vector6d &t_Lambda,
                                          const ::LieAlgebra::Vector6d &t_eta,
                                          const ::LieAlgebra::Vector6d &t_dot_eta,
                                          const ::LieAlgebra::Matrix6d &t_ad_Xi,
                                          const double &t_X)const;



    Eigen::VectorXd backwardODEs(const Eigen::VectorXd &t_y,
                                 const ::LieAlgebra::Vector6d &t_Xi,
                                 const ::LieAlgebra::Vector6d &t_dot_Xi,
                                 const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                 const Eigen::MatrixXd &t_BPhi,
                                 const double &t_X)const;



    typedef Eigen::Matrix<double, 37, 1> TangentKinematicState;

    TangentKinematicState m_forward_tangent_kinematics_state_X0 { TangentKinematicState::Zero() };
    TangentKinematicState m_forward_tangent_kinematics_state_X1 { TangentKinematicState::Zero() };


    Eigen::VectorXd m_tangent_dynamics_state_X0 { Eigen::VectorXd::Zero(49+m_ne) };


    TangentKinematicState tangentKinematicsODEs(const TangentKinematicState &t_state,
                                                const ::LieAlgebra::Vector6d &t_Xi,
                                                const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                                const ::LieAlgebra::Vector6d &t_Delta_Xi,
                                                const ::LieAlgebra::Vector6d &t_Delta_dot_Xi,
                                                const ::LieAlgebra::Vector6d &t_Delta_ddot_Xi) const;


    Eigen::VectorXd tangentDynamicsODEs(const Eigen::VectorXd &t_state,
                                        const ::LieAlgebra::Vector6d &t_Xi,
                                        const ::LieAlgebra::Vector6d &t_dot_Xi,
                                        const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                        const ::LieAlgebra::Vector6d &t_Delta_Xi,
                                        const ::LieAlgebra::Vector6d &t_Delta_dot_Xi,
                                        const ::LieAlgebra::Vector6d &t_Delta_ddot_Xi,
                                        const Eigen::MatrixXd &t_BPhi,
                                        const double &t_X)const;

#ifdef TESTING
public:
    Eigen::MatrixXd m_forward_kinematics_stack {
        Eigen::MatrixXd::Zero(ForwardKinematicState::RowsAtCompileTime, m_number_of_Chebyshev_points)
    };



    Eigen::MatrixXd m_backward_Lambda_Qa_stack {
        Eigen::MatrixXd::Zero(6 + m_ne, m_number_of_Chebyshev_points)
    };

    std::vector<double> m_Chebyshev_points { ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points) };


#endif

};

}   //  namespace CROSP::numerical_integrators::runge_kutta

#endif // RUNGEKUTTAINTEGRATOR_HPP
