#include "CROSP/numerical_integrators/runge_kutta/runge_kutta_integrator.hpp"

#include <boost/numeric/odeint.hpp>


namespace CROSP::numerical_integrators::runge_kutta {


typedef boost::numeric::odeint::runge_kutta_dopri5< Eigen::VectorXd, double,
                                                    Eigen::VectorXd, double,
                                                    boost::numeric::odeint::vector_space_algebra> eigen_stepper;


Eigen::Matrix4d getA(const Eigen::Vector3d t_k){
    Eigen::Matrix4d A;
    A   <<     0  , -t_k(0),  -t_k(1),  -t_k(2),
            t_k(0),     0  ,   t_k(2),  -t_k(1),
            t_k(1), -t_k(2),      0  ,   t_k(0),
            t_k(2),  t_k(1),  -t_k(0),      0  ;

    return A;
}



RungeKuttaIntegrator::RungeKuttaIntegrator(const strain_parameterisation::StrainParameterisation t_strain_parameterisation,
                                           const strain_parameterisation::StrainParameterisation t_strain_parameterisation_Delta,
                                           const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                           std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                           const unsigned int t_number_of_Chebyshev_points)
    :
      m_strain_parameterisation(t_strain_parameterisation),
      m_strain_parameterisation_Delta(t_strain_parameterisation_Delta),
      m_polynomial_representation(t_polynomial_representation),
      m_rod_properties(t_rod_properties),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


std::string RungeKuttaIntegrator::printIntegratorProperties()const
{
    std::stringstream integrator_properties;
    integrator_properties << "Runge Kutta integrator\n" "Adaptive step with initial value of : " << m_dX << " m";

    return integrator_properties.str();

}


void RungeKuttaIntegrator::updateParameterisation(const Eigen::VectorXd &t_qe,
                                                  const Eigen::VectorXd &t_dot_qe,
                                                  const Eigen::VectorXd &t_ddot_qe)
{
    m_qe      = t_qe;
    m_dot_qe  = t_dot_qe;
    m_ddot_qe = t_ddot_qe;
}


void RungeKuttaIntegrator::forwardKinematics()
{
    forwardKinematics(Eigen::Vector4d(1, 0, 0, 0),
                      Eigen::Vector3d::Zero(),
                      Eigen::Vector3d::Zero(),
                      Eigen::Vector3d::Zero(),
                      Eigen::Vector3d::Zero(),
                      Eigen::Vector3d::Zero());
}



void RungeKuttaIntegrator::forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                             const Eigen::Vector3d &t_initial_position,
                                             const Eigen::Vector3d &t_initial_angular_velocity,
                                             const Eigen::Vector3d &t_initial_linear_velocity,
                                             const Eigen::Vector3d &t_initial_angular_acceleration,
                                             const Eigen::Vector3d &t_initial_linear_acceleration)
{

    Eigen::VectorXd forward_integration_state { Eigen::VectorXd::Zero(19) };

    forward_integration_state <<
         t_initial_quaternion,
         t_initial_position,
         t_initial_angular_velocity,
         t_initial_linear_velocity,
         t_initial_angular_acceleration,
         t_initial_linear_acceleration;

    //  Save a track of the state
    m_forward_integration_state_X0 = forward_integration_state;



//  Forward Integration
#ifndef TESTING
    boost::numeric::odeint::integrate_adaptive(
#else
    boost::numeric::odeint::integrate_times(
#endif
                                                eigen_stepper(),
                                               [this](const Eigen::VectorXd &t_y,
                                                      Eigen::VectorXd &t_dyds,
                                                      const double t_X){
        /*  Preprocessing    */
        Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

        //  Get the strains for the rod
        const ::LieAlgebra::Vector6d Xi      = BPhi*m_qe + m_strain_parameterisation.m_constant_strain;
        const ::LieAlgebra::Vector6d dot_Xi  = BPhi*m_dot_qe;
        const ::LieAlgebra::Vector6d ddot_Xi = BPhi*m_ddot_qe;

        t_dyds = this->forwardODEs(t_y, Xi, dot_Xi, ddot_Xi);

                                                },
#ifndef TESTING
                                                forward_integration_state,
                                                0.0,
                                                1.0,
                                                m_dX);
#else
                                                forward_integration_state,
                                                m_Chebyshev_points_reversed.begin(),
                                                m_Chebyshev_points_reversed.end(),
                                                m_dX,
    [&, this, col=0](const Eigen::VectorXd& t_y, const double)mutable{
        m_forward_kinematics_stack.block<ForwardKinematicState::RowsAtCompileTime, 1>(0, col++) = t_y;
    });
#endif

    //  Save a track of the state
    m_forward_integration_state_X1 = forward_integration_state;
}



void RungeKuttaIntegrator::forwardTangentKinematics()
{
    forwardTangentKinematics(Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero());
}



void RungeKuttaIntegrator::updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                                       const Eigen::VectorXd &t_Delta_dot_qe,
                                                       const Eigen::VectorXd &t_Delta_ddot_qe)
{
    m_Delta_qe      = t_Delta_qe;
    m_Delta_dot_qe  = t_Delta_dot_qe;
    m_Delta_ddot_qe = t_Delta_ddot_qe;
}


void RungeKuttaIntegrator::forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                                    const Eigen::Vector3d &t_initial_Delta_position,
                                                    const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                                    const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                                    const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                                    const Eigen::Vector3d &t_initial_Delta_linear_acceleration)
{
    TangentKinematicState tangent_kinematics_state;

    tangent_kinematics_state <<
         m_forward_integration_state_X0,
         t_initial_Delta_orientation,
         t_initial_Delta_position,
         t_initial_Delta_angular_velocity,
         t_initial_Delta_linear_velocity,
         t_initial_Delta_angular_acceleration,
         t_initial_Delta_linear_acceleration;

    //  Save a track of the state
    m_forward_tangent_kinematics_state_X0 = tangent_kinematics_state;



//  Forward Integration
    boost::numeric::odeint::integrate_adaptive(eigen_stepper(),
                                               [this](const Eigen::VectorXd &t_y,
                                                      Eigen::VectorXd &t_dyds,
                                                      const double t_X){
        /*  Preprocessing    */
        Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

        //  Get the strains for the rod
        const ::LieAlgebra::Vector6d Xi      = BPhi*m_qe + m_strain_parameterisation.m_constant_strain;
        const ::LieAlgebra::Vector6d dot_Xi  = BPhi*m_dot_qe;
        const ::LieAlgebra::Vector6d ddot_Xi = BPhi*m_ddot_qe;


        //  Get the strains for the rod
        const ::LieAlgebra::Vector6d Delta_Xi      = BPhi*m_Delta_qe;
        const ::LieAlgebra::Vector6d Delta_dot_Xi  = BPhi*m_Delta_dot_qe;
        const ::LieAlgebra::Vector6d Delta_ddot_Xi = BPhi*m_Delta_ddot_qe;


        t_dyds = this->tangentKinematicsODEs(t_y, Xi, dot_Xi, ddot_Xi, Delta_Xi, Delta_dot_Xi, Delta_ddot_Xi);

                                                },
                                                tangent_kinematics_state,
                                                0.0,
                                                1.0,
                                                m_dX);
    //  Save a track of the state
    m_forward_tangent_kinematics_state_X1 = tangent_kinematics_state;
}


::LieAlgebra::Kinematics RungeKuttaIntegrator::getKinematicsAtTip()const
{
    //  Get the quaternion at the rod tip
    const Eigen::Quaterniond Q(m_forward_integration_state_X1[0],
                               m_forward_integration_state_X1[1],
                               m_forward_integration_state_X1[2],
                               m_forward_integration_state_X1[3]);

    //  Get the quaternion at the rod tip
    const Eigen::Vector3d r(m_forward_integration_state_X1[4],
                            m_forward_integration_state_X1[5],
                            m_forward_integration_state_X1[6]);


    ::LieAlgebra::Kinematics rod_tip_kinematics;


    rod_tip_kinematics.m_pose =
            ::LieAlgebra::SE3Pose( Q, r );


    rod_tip_kinematics.m_twist <<
            m_forward_integration_state_X1.block<6, 1>(7, 0);

    rod_tip_kinematics.m_accelerations <<
            m_forward_integration_state_X1.block<6, 1>(13, 0);

    return rod_tip_kinematics;
}


::LieAlgebra::TangentKinematics RungeKuttaIntegrator::getTangentKinematicsAtTip()const
{
    ::LieAlgebra::TangentKinematics rod_tip_tangent_kinematics;


    rod_tip_tangent_kinematics.m_Delta_zeta =
            m_forward_tangent_kinematics_state_X1.block<6,1>(19,0);

    rod_tip_tangent_kinematics.m_Delta_twist =
            m_forward_tangent_kinematics_state_X1.block<6,1>(25,0);

    rod_tip_tangent_kinematics.m_Delta_acceleration =
            m_forward_tangent_kinematics_state_X1.block<6,1>(31,0);

    return rod_tip_tangent_kinematics;
}






void RungeKuttaIntegrator::backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1)
{
    /*  The backward state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     *  | η |   Ω1, Ω2, Ω3, V1, V2, V3      7-12
     *  | η̇ |   Ω1, Ω2, Ω3, V1, V2, V3     13-18
     *  | Λ |   C1, C2, C3, N1, N2, N3     19-24
     *  | Qa|                              25-(25+ne*na)
     *
     *
     */

    Eigen::VectorXd Q1 = Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension(), 1);

    m_backward_integration_state_X1 <<
                  m_forward_integration_state_X1,
                  t_Lambda_X1,
                  Q1;


    m_backward_integration_state_X0 = m_backward_integration_state_X1;



//  Backward Integration
#ifndef TESTING
    boost::numeric::odeint::integrate_adaptive(
#else
    boost::numeric::odeint::integrate_times(
#endif
                                                eigen_stepper(),
                                               [this](const Eigen::VectorXd& t_y,
                                                      Eigen::VectorXd& t_dyds,
                                                      const double t_X){
        /*  Preprocessing    */
        Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

        //  Get the strains for the rod
        const Eigen::VectorXd Xi      = BPhi*m_qe + m_strain_parameterisation.m_constant_strain;
        const Eigen::VectorXd dot_Xi  = BPhi*m_dot_qe;
        const Eigen::VectorXd ddot_Xi = BPhi*m_ddot_qe;

        t_dyds = this->backwardODEs(t_y, Xi, dot_Xi, ddot_Xi, BPhi, t_X);

                                                    },
#ifndef TESTING
                                                    m_backward_integration_state_X0,
                                                    1.0,
                                                    0.0,
                                                    -m_dX);
#else
                                                    m_backward_integration_state_X0,
                                                    m_Chebyshev_points.begin(),
                                                    m_Chebyshev_points.end(),
                                                    -m_dX,
    [&, this, col=0](const Eigen::VectorXd& t_y, const double)mutable{
        m_backward_Lambda_Qa_stack.block(0, m_number_of_Chebyshev_points-1-(col++), 6+m_ne, 1) = t_y.block(19, 0, 6+m_ne, 1);
    });
#endif



}





void RungeKuttaIntegrator::backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
{

    /*  The state has the form
     *  | Q  |   w, x, y, z                       0-3
     *  | r  |   x, y, z                          4-6
     *  | η  |   Ω1, Ω2, Ω3, V1, V2, V3          7-12
     *  | η̇  |   Ω1, Ω2, Ω3, V1, V2, V3         13-18
     *  | ∆ζ |  ∆K1, ∆K2, ∆K3, ∆Γ1, ∆Γ2, ∆Γ3    19-24
     *  | ∆η |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   25-30
     *  | ∆η̇ |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   31-36
     *  | Λ  |   C1, C2, C3, N1, N2, N3          37-42
     *  | ∆Λ |   C1, C2, C3, N1, N2, N3          44+ne-49+ne
     *  | ∆Qa|                                   50+ne-50+2ne
     */

    Eigen::VectorXd Delta_Q1 = Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension(), 1);

    Eigen::VectorXd Lambda_X1 = m_backward_integration_state_X1.block<6, 1>(19,0);

    m_tangent_dynamics_state_X0 <<
                  m_forward_tangent_kinematics_state_X1,
                  Lambda_X1,
                  t_Delta_Lambda_X1,
                  Delta_Q1;



//  Backward Integration
    boost::numeric::odeint::integrate_adaptive(eigen_stepper(),
                                               [this](const Eigen::VectorXd& t_y,
                                                      Eigen::VectorXd& t_dyds,
                                                      const double t_X){
        /*  Preprocessing    */
        Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

        //  Get the strains for the rod
        const ::LieAlgebra::Vector6d Xi      = BPhi*m_qe + m_strain_parameterisation.m_constant_strain;
        const ::LieAlgebra::Vector6d dot_Xi  = BPhi*m_dot_qe;
        const ::LieAlgebra::Vector6d ddot_Xi = BPhi*m_ddot_qe;


        //  Get the strains for the rod
        const ::LieAlgebra::Vector6d Delta_Xi      = BPhi*m_Delta_qe;
        const ::LieAlgebra::Vector6d Delta_dot_Xi  = BPhi*m_Delta_dot_qe;
        const ::LieAlgebra::Vector6d Delta_ddot_Xi = BPhi*m_Delta_ddot_qe;

        t_dyds = this->tangentDynamicsODEs(t_y,
                                           Xi, dot_Xi, ddot_Xi,
                                           Delta_Xi, Delta_dot_Xi, Delta_ddot_Xi,
                                           BPhi, t_X);

                                                      },
                                                      m_tangent_dynamics_state_X0,
                                                      1.0,
                                                      0.0,
                                                      -m_dX);

}



::LieAlgebra::Vector6d RungeKuttaIntegrator::getLambdaAtBase()const
{
    return m_backward_integration_state_X0.block<6, 1>(19, 0);
}



LieAlgebra::Vector6d RungeKuttaIntegrator::getDeltaLambdaAtBase()const
{
    return m_tangent_dynamics_state_X0.block<6, 1>(43,0);
}


::LieAlgebra::Vector6d RungeKuttaIntegrator::getQaAtBase()const
{
    return m_backward_integration_state_X0.block(25, 0, m_ne, 1);
}



LieAlgebra::Vector6d RungeKuttaIntegrator::getDeltaQaAtBase()const
{
    return m_tangent_dynamics_state_X0.block(49, 0, m_ne, 1);
}


void RungeKuttaIntegrator::updateIntegrationDomain(const double &t_rod_lenght)
{
    m_rod_length = t_rod_lenght;
}






Eigen::MatrixXd RungeKuttaIntegrator::getRodPositions() const
{

    /*  The forward state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     */
    Eigen::VectorXd y = m_forward_integration_state_X0.block<7, 1>(0, 0);


    Eigen::MatrixXd rod_shapes_stack { Eigen::MatrixXd::Zero(3, m_number_of_Chebyshev_points) };

//  Forward Integration with observer
    boost::numeric::odeint::integrate_times(eigen_stepper(),
                                            [this](const Eigen::VectorXd& t_y, Eigen::VectorXd& t_dydX, const double t_s){

        /*  Preprocessing    */
        Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_s);

        //  Get the strains for the rod
        const Eigen::VectorXd Xi = BPhi*m_qe + m_strain_parameterisation.m_constant_strain;


        t_dydX = this->forwardStaticODEs(t_y, Xi);


    },
    y,
    m_Chebyshev_points_reversed.begin(),
    m_Chebyshev_points_reversed.end(),
    m_dX,
    [&, col=0](const Eigen::VectorXd& t_y, const double)mutable{
        rod_shapes_stack.block<3, 1>(0, col++) = t_y.block<3,1>(4, 0);
    });


    return rod_shapes_stack;
}






RungeKuttaIntegrator::PoseState RungeKuttaIntegrator::forwardStaticODEs(const RungeKuttaIntegrator::PoseState &t_state,
                                                                        const ::LieAlgebra::Vector6d &t_Xi) const
{
    //  Decompose the strain
    const Eigen::Vector3d k = t_Xi.block<3,1>(0,0);
    const Eigen::Vector3d gamma = t_Xi.block<3,1>(3,0);


    /*  The state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     */


    //  Unpack state vector
    const Eigen::Quaterniond Q(t_state[0], t_state[1],t_state[2], t_state[3]);
    [[maybe_unused]] const Eigen::Vector3d r = t_state.block<3,1>(4,0);


    const Eigen::Matrix3d R = Q.toRotationMatrix();


    //  Actual ODE
    const Eigen::Vector4d Q_prime = 0.5*getA(k)*Eigen::Vector4d(Q.w(),Q.x(),Q.y(),Q.z());
    const Eigen::Vector3d r_prime = R*gamma;


    //  Packing state vector derivative
    RungeKuttaIntegrator::PoseState dyds;
    dyds << Q_prime,
            r_prime;

    dyds *= m_rod_length;

    return dyds;

}

RungeKuttaIntegrator::ForwardKinematicState RungeKuttaIntegrator::forwardODEs(const RungeKuttaIntegrator::ForwardKinematicState &t_state,
                                                                              const ::LieAlgebra::Vector6d &t_Xi,
                                                                              const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                                              const ::LieAlgebra::Vector6d &t_ddot_Xi) const
{
    /*  The state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     *  | η |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 7-12
     *  | η̇ |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 13-18
     */


    //  Unpack state vector
    const RungeKuttaIntegrator::PoseState g = t_state.block<7,1>(0,0);
    const ::LieAlgebra::Vector6d eta = t_state.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d eta_dot = t_state.block<6,1>(13,0);


    const auto ad_Xi = ::LieAlgebra::ad(t_Xi);
    const auto ad_dot_Xi = ::LieAlgebra::ad(t_dot_Xi);

    //  Actual ODE
    const RungeKuttaIntegrator::PoseState g_prime = forwardStaticODEs(g, t_Xi);
    const ::LieAlgebra::Vector6d eta_prime = - ad_Xi*eta + t_dot_Xi;
    const ::LieAlgebra::Vector6d eta_dot_prime = - ad_Xi*eta_dot - ad_dot_Xi*eta + t_ddot_Xi;



    RungeKuttaIntegrator::ForwardKinematicState dydx;

    //  Packing state vector derivative
    dydx << g_prime,
            eta_prime,
            eta_dot_prime;

    dydx *= m_rod_length;

    return dydx;
}






::LieAlgebra::Vector6d RungeKuttaIntegrator::getLambdaPrime(const Eigen::Quaterniond &t_Q,
                                                            const ::LieAlgebra::Vector6d &t_Lambda,
                                                            const ::LieAlgebra::Vector6d &t_eta,
                                                            const ::LieAlgebra::Vector6d &t_dot_eta,
                                                            const ::LieAlgebra::Matrix6d &t_ad_Xi,
                                                            const double &t_X)const
{
    if(t_X == 0 or t_X == 1.0)
        return ::LieAlgebra::Vector6d::Zero();


    //  Some needed variables
    const Eigen::Matrix3d R = t_Q.toRotationMatrix();
    ::LieAlgebra::Vector6d F_bar = ::LieAlgebra::Vector6d::Zero();
    F_bar.block<3, 1>(3, 0) = R.transpose() * m_rod_properties->distributedGravitationalForce();


    const auto M = m_rod_properties->m_M;

    const Eigen::VectorXd Lambda_prime = t_ad_Xi.transpose()*t_Lambda
                                            + M*t_dot_eta
                                            - ::LieAlgebra::ad(t_eta).transpose()*M*t_eta
                                            - F_bar;

    return Lambda_prime;
}



Eigen::VectorXd RungeKuttaIntegrator::backwardODEs(const Eigen::VectorXd &t_y,
                                                   const ::LieAlgebra::Vector6d &t_Xi,
                                                   const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                   const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                                   const Eigen::MatrixXd &t_BPhi,
                                                   const double &t_X)const
{
    /*  The state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     *  | η |   Ω1, Ω2, Ω3, V1, V2, V3      7-12
     *  | η̇ |   Ω1, Ω2, Ω3, V1, V2, V3     13-18
     *  | Λ |   C1, C2, C3, N1, N2, N3     19-24
     *  | Qa|                              25-(25+ne*na)
     */

    //  Unpack state vector
    const ForwardKinematicState kinematic_state = t_y.block<19,1>(0,0);

    const ::LieAlgebra::Vector6d eta = t_y.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d dot_eta = t_y.block<6,1>(13,0);
    const ::LieAlgebra::Vector6d Lambda = t_y.block<6,1>(19,0);




    const auto ad_Xi = ::LieAlgebra::ad(t_Xi);



    //  Actual ODE
    const ForwardKinematicState kinematic_state_prime = forwardODEs(kinematic_state, t_Xi, t_dot_Xi, t_ddot_Xi);


    const Eigen::VectorXd Lambda_prime = getLambdaPrime(Eigen::Quaterniond(t_y[0], t_y[1],t_y[2], t_y[3]),
                                                                           Lambda, eta, dot_eta, ad_Xi, t_X);
    const Eigen::VectorXd Qa_prime = - t_BPhi.transpose()*Lambda;

    //  Packing state vector derivative
    Eigen::VectorXd dyds(25 + m_ne);
    dyds <<  kinematic_state_prime,
             Lambda_prime,
             Qa_prime;

    dyds *= m_rod_length;

    return dyds;
}



RungeKuttaIntegrator::TangentKinematicState RungeKuttaIntegrator::tangentKinematicsODEs(const TangentKinematicState &t_state,
                                                                                        const ::LieAlgebra::Vector6d &t_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_Delta_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_Delta_dot_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_Delta_ddot_Xi) const
{
    /*  The state has the form
     *  | Q |   w, x, y, z                      0-3
     *  | r |   x, y, z                         4-6
     *  | η |   Ω1, Ω2, Ω3, V1, V2, V3         7-12
     *  | η̇ |   Ω1, Ω2, Ω3, V1, V2, V3         13-18
     *  | ∆ζ |  ∆K1, ∆K2, ∆K3, ∆Γ1, ∆Γ2, ∆Γ3    19-24
     *  | ∆η |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   25-30
     *  | ∆η̇ |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   31-37
     */

    //  Unpack state vector
    const ForwardKinematicState kinematic_state = t_state.block<19,1>(0,0);
    const ::LieAlgebra::Vector6d eta = t_state.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d dot_eta = t_state.block<6,1>(13,0);

    const ::LieAlgebra::Vector6d Delta_zeta    = t_state.block<6,1>(19,0);
    const ::LieAlgebra::Vector6d Delta_eta     = t_state.block<6,1>(25,0);
    const ::LieAlgebra::Vector6d Delta_dot_eta = t_state.block<6,1>(31,0);



    const auto ad_Xi = ::LieAlgebra::ad(t_Xi);

    const auto ad_Delta_dot_Xi = ::LieAlgebra::ad(t_Delta_dot_Xi);

    const auto ad_eta = ::LieAlgebra::ad(eta);
    const auto ad_dot_eta = ::LieAlgebra::ad(dot_eta);


    const ForwardKinematicState kinematic_state_prime = forwardODEs(kinematic_state, t_Xi, t_dot_Xi, t_ddot_Xi);
    const ::LieAlgebra::Vector6d Delta_zeta_prime =
            - ad_Xi*Delta_zeta + t_Delta_Xi;
    const ::LieAlgebra::Vector6d Delta_eta_prime =
            - ad_Xi*Delta_eta + ad_eta*t_Delta_Xi + t_Delta_dot_Xi;
    const ::LieAlgebra::Vector6d Delta_dot_eta_prime =
            - ad_Xi*Delta_dot_eta - ad_Delta_dot_Xi*Delta_eta + ad_eta*t_Delta_dot_Xi + ad_dot_eta*t_Delta_Xi + t_Delta_ddot_Xi;


    //  Packing state vector derivative
    TangentKinematicState dydx;
    dydx <<  kinematic_state_prime,
             Delta_zeta_prime,
             Delta_eta_prime,
             Delta_dot_eta_prime;

    dydx *= m_rod_length;

    return dydx;
}


Eigen::VectorXd RungeKuttaIntegrator::tangentDynamicsODEs(const Eigen::VectorXd &t_state,
                                                          const ::LieAlgebra::Vector6d &t_Xi,
                                                          const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                          const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                                          const ::LieAlgebra::Vector6d &t_Delta_Xi,
                                                          const ::LieAlgebra::Vector6d &t_Delta_dot_Xi,
                                                          const ::LieAlgebra::Vector6d &t_Delta_ddot_Xi,
                                                          const Eigen::MatrixXd &t_BPhi,
                                                          const double &t_X)const
{

    /*  The state has the form
     *  | Q  |   w, x, y, z                       0-3
     *  | r  |   x, y, z                          4-6
     *  | η  |   Ω1, Ω2, Ω3, V1, V2, V3          7-12
     *  | η̇  |   Ω1, Ω2, Ω3, V1, V2, V3         13-18
     *  | ∆ζ |  ∆K1, ∆K2, ∆K3, ∆Γ1, ∆Γ2, ∆Γ3    19-24
     *  | ∆η |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   25-30
     *  | ∆η̇ |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   31-36
     *  | Λ  |   C1, C2, C3, N1, N2, N3          37-42
     *  | ∆Λ |   C1, C2, C3, N1, N2, N3          43-48
     *  | ∆Qa|                                   49-49+ne
     */


    //  Unpack state vector
    const TangentKinematicState tangent_kinematic_state = t_state.block<37,1>(0,0);
    const ::LieAlgebra::Vector6d eta = t_state.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d dot_eta = t_state.block<6,1>(13,0);

    const ::LieAlgebra::Vector6d Delta_zeta    = t_state.block<6,1>(19,0);
    const ::LieAlgebra::Vector6d Delta_eta     = t_state.block<6,1>(25,0);
    const ::LieAlgebra::Vector6d Delta_dot_eta = t_state.block<6,1>(31,0);


    const ::LieAlgebra::Vector6d Lambda = t_state.block<6, 1>(37,0);
    const ::LieAlgebra::Vector6d Delta_Lambda = t_state.block<6, 1>(43,0);



    const auto ad_Xi = ::LieAlgebra::ad(t_Xi);
    const auto ad_Delta_Xi = ::LieAlgebra::ad(t_Delta_Xi);


    const auto ad_eta = ::LieAlgebra::ad(eta);
    const auto ad_Delta_eta = ::LieAlgebra::ad(Delta_eta);

    const auto M = m_rod_properties->m_M;








    const Eigen::Matrix3d R = Eigen::Quaterniond(t_state[0], t_state[1],t_state[2], t_state[3]).toRotationMatrix();

    Eigen::Vector3d Delta_rotation = Delta_zeta.block<3,1>(0, 0);

    Eigen::Vector3d Delta_N_bar = ::LieAlgebra::skew( Delta_rotation ).transpose()
                                    *R.transpose()
                                    *m_rod_properties->distributedGravitationalForce();

    ::LieAlgebra::Vector6d Delta_F_bar = ::LieAlgebra::Vector6d::Zero();
    Delta_F_bar.block<3,1>(3, 0) = Delta_N_bar;





    const TangentKinematicState tangent_kinematic_state_prime = tangentKinematicsODEs(tangent_kinematic_state,
                                                                                      t_Xi, t_dot_Xi, t_ddot_Xi,
                                                                                      t_Delta_Xi, t_Delta_dot_Xi, t_Delta_ddot_Xi);
    const Eigen::VectorXd Lambda_prime = getLambdaPrime(Eigen::Quaterniond(t_state[0], t_state[1],t_state[2], t_state[3]),
                                                                           Lambda, eta, dot_eta, ad_Xi, t_X);
    const ::LieAlgebra::Vector6d Delta_Lambda_prime =
            M*Delta_dot_eta - ad_eta.transpose()*M*Delta_eta - ad_Delta_eta.transpose()*M*eta + ad_Xi.transpose()*Delta_Lambda + ad_Delta_Xi.transpose()*Lambda - Delta_F_bar;

    const Eigen::VectorXd Delta_Qa_prime = -t_BPhi.transpose()*Delta_Lambda;


    //  Packing state vector derivative
    Eigen::VectorXd dydx(49 + m_ne);
    dydx <<  tangent_kinematic_state_prime,
             Lambda_prime,
             Delta_Lambda_prime,
             Delta_Qa_prime;

    dydx *= m_rod_length;

    return dydx;
}




}   //  namespace CROSP::numerical_integrators::runge_kutta
