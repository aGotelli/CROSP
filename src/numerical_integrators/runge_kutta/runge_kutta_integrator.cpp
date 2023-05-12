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
                                           const rod_properties::RodProperties t_rod_properties,
                                           const unsigned int t_number_of_Chebyshev_points)
    :
      m_strain_parameterisation(t_strain_parameterisation),
      m_strain_parameterisation_Delta(t_strain_parameterisation_Delta),
      m_polynomial_representation(t_polynomial_representation),
      m_rod_properties(t_rod_properties),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


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
    boost::numeric::odeint::integrate_adaptive(eigen_stepper(),
                                               [this](const Eigen::VectorXd& t_y,
                                                       Eigen::VectorXd& t_dyds,
                                                       const double t_s){
                                                            this->forwardODEs(t_y,
                                                                              t_dyds,
                                                                              t_s);
                                                       },
                                                forward_integration_state,
                                                0.0,
                                                1.0,
                                                m_dX);
    //  Save a track of the state
    m_forward_integration_state_X1 = forward_integration_state;
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

    m_backward_integration_state_X0 <<
                  m_forward_integration_state_X1,
                  t_Lambda_X1,
                  Q1;



//  Backward Integration
    boost::numeric::odeint::integrate_adaptive(eigen_stepper(),
                                               [this](const Eigen::VectorXd& t_y,
                                                      Eigen::VectorXd& t_dyds,
                                                      const double t_s){
                                                        this->backwardODEs(t_y,
                                                                           t_dyds,
                                                                           t_s);
                                                      },
                                                      m_backward_integration_state_X0,
                                                      1.0,
                                                      0.0,
                                                      -m_dX);


}


Eigen::MatrixXd RungeKuttaIntegrator::getRodPositions() const
{

    /*  The forward state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     *  | η |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 7-12
     *  | η̇ |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 13-18
     */
    Eigen::VectorXd y = m_forward_integration_state_X0.block<7, 1>(0, 0);


    Eigen::MatrixXd rod_shapes_stack { Eigen::MatrixXd::Zero(3, m_number_of_Chebyshev_points) };

//  Forward Integration with observer
    boost::numeric::odeint::integrate_times(eigen_stepper(),
                                            [this](const Eigen::VectorXd& t_y, Eigen::VectorXd& t_dyds, const double t_s){

        /*  Preprocessing    */

            Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_s);

        //  Get the strains for the rod
            const Eigen::VectorXd xi = BPhi*m_qe + m_strain_parameterisation.m_constant_strain;


        //  Decompose the strain
            const Eigen::Vector3d k = xi.block<3,1>(0,0);
            const Eigen::Vector3d gamma = xi.block<3,1>(3,0);



        /*  The state has the form
         *  | Q |   w, x, y, z                  0-3
         *  | r |   x, y, z                     4-6
         *  | η |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 7-12
         *  | η̇ |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 13-18
         */


        //  Unpack state vector
            const Eigen::Quaterniond Q(t_y[0], t_y[1],t_y[2], t_y[3]);
            [[maybe_unused]] const Eigen::Vector3d r = t_y.block<3,1>(4,0);


            const Eigen::Matrix3d R = Q.toRotationMatrix();


        //  Actual ODE

            const Eigen::Vector4d Q_prime = 0.5*getA(k)*Eigen::Vector4d(Q.w(),Q.x(),Q.y(),Q.z());
            const Eigen::Vector3d r_prime = R*gamma;

        //  Packing state vector derivative
            t_dyds.resize( t_y.rows() );
            t_dyds << Q_prime,
                      r_prime;

            t_dyds *= m_rod_length;

    },
    y,
    m_Chebyshev_points_reversed.begin(),
    m_Chebyshev_points_reversed.end(),
    m_dX,
    [&, col=0](const Eigen::VectorXd& t_y, const double)mutable{
        std::cout << "col : " << col << "\n";
        rod_shapes_stack.block<3, 1>(0, col++) = t_y.block<3,1>(4, 0);
    });


    return rod_shapes_stack;
}









void RungeKuttaIntegrator::forwardODEs(const Eigen::VectorXd &t_y,
                                       Eigen::VectorXd &t_dyds,
                                       const double t_s) const
{
/*  Preprocessing    */

    Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_s);

//  Get the strains for the rod
    const Eigen::VectorXd xi      = BPhi*m_qe + m_strain_parameterisation.m_constant_strain;
    const Eigen::VectorXd xi_dot  = BPhi*m_dot_qe;
    const Eigen::VectorXd xi_ddot = BPhi*m_ddot_qe;

//  Decompose the strain
    const Eigen::Vector3d k = xi.block<3,1>(0,0);
    const Eigen::Vector3d gamma = xi.block<3,1>(3,0);



/*  The state has the form
 *  | Q |   w, x, y, z                  0-3
 *  | r |   x, y, z                     4-6
 *  | η |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 7-12
 *  | η̇ |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 13-18
 */


//  Unpack state vector
    const Eigen::Quaterniond Q(t_y[0], t_y[1],t_y[2], t_y[3]);
    [[maybe_unused]] const Eigen::Vector3d r = t_y.block<3,1>(4,0);
    const ::LieAlgebra::Vector6d eta = t_y.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d eta_dot = t_y.block<6,1>(13,0);

    const Eigen::Matrix3d R = Q.toRotationMatrix();


//  Actual ODE

    const Eigen::Vector4d Q_prime = 0.5*getA(k)*Eigen::Vector4d(Q.w(),Q.x(),Q.y(),Q.z());
    const Eigen::Vector3d r_prime = R*gamma;
    const ::LieAlgebra::Vector6d eta_prime = - ::LieAlgebra::ad(xi)*eta + xi_dot;
    const ::LieAlgebra::Vector6d eta_dot_prime = - ::LieAlgebra::ad(xi)*eta_dot - ::LieAlgebra::ad(xi_dot)*eta + xi_ddot;

//  Packing state vector derivative
    t_dyds.resize( t_y.rows() );
    t_dyds << Q_prime,
              r_prime,
              eta_prime,
              eta_dot_prime;

    t_dyds *= m_rod_length;

}



void RungeKuttaIntegrator::backwardODEs(const Eigen::VectorXd &t_y,
                                        Eigen::VectorXd &t_dyds,
                                        const double t_s)const
{
    /*  Preprocessing    */
    Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_s);

    //  Get the strains for the rod
    const Eigen::VectorXd xi      = BPhi*m_qe + m_strain_parameterisation.m_constant_strain;
    const Eigen::VectorXd xi_dot  = BPhi*m_dot_qe;
    const Eigen::VectorXd xi_ddot = BPhi*m_ddot_qe;

    //std::cout << "xi : " << xi.transpose() << std::endl;

//  Decompose the strain
    const Eigen::Vector3d k = xi.block<3,1>(0,0);
    const Eigen::Vector3d gamma = xi.block<3,1>(3,0);


/*  The state has the form
 *  | Q |   w, x, y, z                  0-3
 *  | r |   x, y, z                     4-6
 *  | η |   Ω1, Ω2, Ω3, V1, V2, V3      7-12
 *  | η̇ |   Ω1, Ω2, Ω3, V1, V2, V3     13-18
 *  | Λ |   C1, C2, C3, N1, N2, N3     19-24
 *  | Qa|                              25-(25+ne*na)
 */

    unsigned int ne = m_polynomial_representation.getCoordinatesDimension();


    Eigen::MatrixXd m_M = m_rod_properties.m_M;

    //  Unpack state vector
    const Eigen::Quaterniond Q(t_y[0], t_y[1],t_y[2], t_y[3]);
    [[maybe_unused]] const Eigen::Vector3d r = t_y.block<3,1>(4,0);
    const ::LieAlgebra::Vector6d eta = t_y.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d eta_dot = t_y.block<6,1>(13,0);
    const ::LieAlgebra::Vector6d Lambda = t_y.block<6,1>(19,0);
    const Eigen::VectorXd Qa = t_y.block(25,0, ne,1);

    //std::cout << "Lambda : " << Lambda.transpose() << std::endl;

    const Eigen::Matrix3d R = Q.toRotationMatrix();
    const Eigen::VectorXd F_bar = Eigen::VectorXd::Zero(6, 1);



//  Actual ODE
    const Eigen::Vector4d Q_prime = 0.5*getA(k)*Eigen::Vector4d(Q.w(),Q.x(),Q.y(),Q.z());
    const Eigen::Vector3d r_prime = R*gamma;
    const ::LieAlgebra::Vector6d eta_prime = -::LieAlgebra::ad(xi)*eta + xi_dot;
    const ::LieAlgebra::Vector6d eta_dot_prime = -::LieAlgebra::ad(xi)*eta_dot - ::LieAlgebra::ad(xi_dot)*eta + xi_ddot;
    const Eigen::VectorXd Lambda_prime = ::LieAlgebra::ad(xi).transpose()*Lambda + m_M*eta_dot - ::LieAlgebra::ad(eta).transpose()*m_M*eta - F_bar;
    const Eigen::VectorXd Qa_prime = -BPhi.transpose()*Lambda;

//  Packing state vector derivative
    t_dyds.resize( t_y.rows() );
    t_dyds  <<  Q_prime,
                r_prime,
                eta_prime,
                eta_dot_prime,
                Lambda_prime,
                Qa_prime;

    t_dyds *= m_rod_length;
}


}   //  namespace CROSP::numerical_integrators::runge_kutta
