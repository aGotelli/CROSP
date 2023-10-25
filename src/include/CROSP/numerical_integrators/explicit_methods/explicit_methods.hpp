#ifndef EXPLICIT_METHOTDS_HPP
#define EXPLICIT_METHOTDS_HPP



//#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include <boost/numeric/odeint.hpp>


#include "ATORS/ATORS.hpp"






namespace CROSP::numerical_integrators::explicit_methods {



typedef boost::numeric::odeint::runge_kutta_dopri5< Eigen::VectorXd, double,
                                                    Eigen::VectorXd, double,
                                                    boost::numeric::odeint::vector_space_algebra> runge_kutta_dopri5;

typedef boost::numeric::odeint::runge_kutta4< Eigen::VectorXd, double,
                                                    Eigen::VectorXd, double,
                                                    boost::numeric::odeint::vector_space_algebra> runge_kutta4;


typedef boost::numeric::odeint::bulirsch_stoer< Eigen::VectorXd, double,
                                                    Eigen::VectorXd, double,
                                                    boost::numeric::odeint::vector_space_algebra> bulirsch_stoer;

typedef boost::numeric::odeint::runge_kutta_cash_karp54< Eigen::VectorXd, double,
                                                    Eigen::VectorXd, double,
                                                    boost::numeric::odeint::vector_space_algebra> runge_kutta_cash_karp54;

typedef boost::numeric::odeint::euler< Eigen::VectorXd, double,
                                                    Eigen::VectorXd, double,
                                                    boost::numeric::odeint::vector_space_algebra> euler;


typedef boost::numeric::odeint::modified_midpoint< Eigen::VectorXd, double,
                                                    Eigen::VectorXd, double,
                                                    boost::numeric::odeint::vector_space_algebra> modified_midpoint;



std::string extractBoostSignatureFromtypeid(std::string t_method_type);

/*!
 * \brief The RungeKuttaIntegrator class provides a numeric integrator for the Cosserat IDM and TIDM using the standard Runge-Kutta integrator
 */
class ExplicitIntegrationODEs
{
public:
    ExplicitIntegrationODEs()=default;

    ExplicitIntegrationODEs(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                       std::shared_ptr<const rod_properties::RodProperties> t_rod_properties);

    //  Define a compiled time known dimension for the state of the pose (quaternion + position)
    typedef Eigen::Matrix<double, 7, 1> PoseState;

    typedef Eigen::Matrix<double, 19, 1> ForwardKinematicState;

    typedef Eigen::Matrix<double, 37, 1> TangentKinematicState;





    void forwardStaticODEs(const Eigen::VectorXd &t_y,
                           Eigen::VectorXd &t_dyds,
                           const double t_X) const;




    void forwardODEs(const Eigen::VectorXd &t_y,
                           Eigen::VectorXd &t_dyds,
                           const double t_X) const;



    void backwardODEs(const Eigen::VectorXd &t_y,
                                              Eigen::VectorXd &t_dyds,
                                              const double t_X)const;

    [[deprecated("This function does not work.. there is a bug in Delta_dot_V")]]
    void tangentKinematicsODEs(const Eigen::VectorXd &t_y,
                                              Eigen::VectorXd &t_dyds,
                                              const double t_X)const;


    void tangentDynamicsODEs(const Eigen::VectorXd &t_y,
                                              Eigen::VectorXd &t_dyds,
                                              const double t_X)const;





    void distributedActuationODE(const Eigen::VectorXd &t_y,
                                 Eigen::VectorXd &t_dyds,
                                 const double t_X)const;


    void PhiTPhiODE(const Eigen::VectorXd &t_y,
                    Eigen::VectorXd &t_dyds,
                    const double t_X)const;


private:

    PoseState forwardStaticStep(const PoseState &t_state,
                                const ::LieAlgebra::Vector6d &t_Xi) const;


    ForwardKinematicState forwardStep(const ForwardKinematicState &t_state,
                                      const ::LieAlgebra::Vector6d &t_Xi,
                                      const ::LieAlgebra::Vector6d &t_dot_Xi,
                                      const ::LieAlgebra::Vector6d &t_ddot_Xi) const;



    ::LieAlgebra::Vector6d getLambdaPrime(const Eigen::Quaterniond &t_Q,
                                          const ::LieAlgebra::Vector6d &t_Lambda,
                                          const ::LieAlgebra::Vector6d &t_eta,
                                          const ::LieAlgebra::Vector6d &t_dot_eta,
                                          const ::LieAlgebra::Matrix6d &t_ad_Xi)const;





    Eigen::VectorXd backwardStep(const Eigen::VectorXd &t_state,
                                 const ::LieAlgebra::Vector6d &t_Xi,
                                 const ::LieAlgebra::Vector6d &t_dot_Xi,
                                 const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                 const Eigen::MatrixXd &t_BPhi)const;






    TangentKinematicState tangentKinematicsStep(const TangentKinematicState &t_state,
                                                const ::LieAlgebra::Vector6d &t_Xi,
                                                const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                                const ::LieAlgebra::Vector6d &t_Delta_Xi,
                                                const ::LieAlgebra::Vector6d &t_Delta_dot_Xi,
                                                const ::LieAlgebra::Vector6d &t_Delta_ddot_Xi) const;



    Eigen::VectorXd tangentDynamicsStep(const Eigen::VectorXd &t_state,
                                        const ::LieAlgebra::Vector6d &t_Xi,
                                        const ::LieAlgebra::Vector6d &t_dot_Xi,
                                        const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                        const ::LieAlgebra::Vector6d &t_Delta_Xi,
                                        const ::LieAlgebra::Vector6d &t_Delta_dot_Xi,
                                        const ::LieAlgebra::Vector6d &t_Delta_ddot_Xi,
                                        const Eigen::MatrixXd &t_BPhi)const;






    std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;


    polynomial_representation::PolynomialRepresentation m_polynomial_representation;
    ::LieAlgebra::Vector6d m_constant_strain { ::LieAlgebra::Vector6d::Unit(3) };



protected:



    double m_rod_length { m_rod_properties->m_rod_dimensions.m_L };

    unsigned int m_generalised_coordinates_dimension { m_polynomial_representation.getCoordinatesDimension() };


    Eigen::VectorXd m_qe { Eigen::VectorXd::Zero(m_generalised_coordinates_dimension) };
    Eigen::VectorXd m_dot_qe { Eigen::VectorXd::Zero(m_generalised_coordinates_dimension) };
    Eigen::VectorXd m_ddot_qe { Eigen::VectorXd::Zero(m_generalised_coordinates_dimension) };


    Eigen::VectorXd m_Delta_qe { Eigen::VectorXd::Zero(m_generalised_coordinates_dimension) };
    Eigen::VectorXd m_Delta_dot_qe { Eigen::VectorXd::Zero(m_generalised_coordinates_dimension) };
    Eigen::VectorXd m_Delta_ddot_qe { Eigen::VectorXd::Zero(m_generalised_coordinates_dimension) };






    Eigen::VectorXd m_forward_integration_state_X0 { Eigen::VectorXd::Zero(19) };
    Eigen::VectorXd m_forward_integration_state_X1 { Eigen::VectorXd::Zero(19) };



    unsigned int m_backward_state_size { 25+m_generalised_coordinates_dimension };
    Eigen::VectorXd m_backward_integration_state_X0 { Eigen::VectorXd::Zero(25+m_generalised_coordinates_dimension) };
    Eigen::VectorXd m_backward_integration_state_X1 { Eigen::VectorXd::Zero(25+m_generalised_coordinates_dimension) };


    Eigen::VectorXd m_forward_tangent_kinematics_state_X0 { Eigen::VectorXd::Zero(37) };
    Eigen::VectorXd m_forward_tangent_kinematics_state_X1 { Eigen::VectorXd::Zero(37) };


    Eigen::VectorXd m_tangent_dynamics_state_X0 { Eigen::VectorXd::Zero(49+m_generalised_coordinates_dimension) };
    Eigen::VectorXd m_tangent_dynamics_state_X1 { Eigen::VectorXd::Zero(49+m_generalised_coordinates_dimension) };


    //::ATORS::tendon_driven_actuation::TendonDrivenActuation m_distributed_actuation;
    ::ATORS::distributed_actuation::DistributedActuationUptr m_distributed_actuation;

    Eigen::VectorXd m_Qad_X0 { Eigen::VectorXd::Zero(m_generalised_coordinates_dimension) };


    std::function<void()> m_integrate_Qad { [](){} };
};









template<class Method=runge_kutta_dopri5>
class ExplicitIntegrator : public ExplicitIntegrationODEs {

public:

    ExplicitIntegrator(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                       const unsigned int t_number_of_Chebyshev_points,
                       std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
        : ExplicitIntegrationODEs(t_polynomial_representation,
                                  t_rod_properties),
          m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
    {}


    void addInternalActuation(::ATORS::distributed_actuation::DistributedActuationUptr t_distributed_actuation)
    {
        m_distributed_actuation = std::move(t_distributed_actuation);

        m_integrate_Qad = [this]() {
            boost::numeric::odeint::integrate_adaptive(stepper,
                                                       [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                        {this->distributedActuationODE(t_y, t_dyds, t_s);},
                                                         m_Qad_X0,
                                                         1.0,
                                                         0.0,
                                                         -m_dX);
        };
    }





    std::string printIntegratorProperties()const
    {

        std::string method_type =  extractBoostSignatureFromtypeid( typeid(stepper).name() );


        std::stringstream integrator_properties;
        integrator_properties << "Stepper integrator\n"
                                 "      Method : " << method_type << "   Order : " << stepper.order() << "\n"
                                 "      Adaptive step with initial value of : " << m_dX << "\n"
                                 "Integration domain : [0, 1]\n"
                                 "Observerving using : " << m_number_of_Chebyshev_points << " Chebyshev points";

        return integrator_properties.str();

    }


    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                                      const Eigen::VectorXd &t_dot_qe,
                                                      const Eigen::VectorXd &t_ddot_qe)
    {
        m_qe      = t_qe;
        m_dot_qe  = t_dot_qe;
        m_ddot_qe = t_ddot_qe;
    }



    void forwardKinematics()
    {
        forwardKinematics(Eigen::Vector4d(1, 0, 0, 0),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero());
    }


    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                           const Eigen::Vector3d &t_initial_position,
                           const Eigen::Vector3d &t_initial_angular_velocity,
                           const Eigen::Vector3d &t_initial_linear_velocity,
                           const Eigen::Vector3d &t_initial_angular_acceleration,
                           const Eigen::Vector3d &t_initial_linear_acceleration)
    {

        ExplicitIntegrationODEs::ForwardKinematicState forward_integration_state;
        forward_integration_state << t_initial_quaternion,
                t_initial_position,
                t_initial_angular_velocity,
                t_initial_linear_velocity,
                t_initial_angular_acceleration,
                t_initial_linear_acceleration;

        //  Save a track of the state
        m_forward_integration_state_X0 = forward_integration_state;


        boost::numeric::odeint::integrate_adaptive(stepper,
                                                   [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                        {this->forwardODEs(t_y, t_dyds, t_s);},
                                                   forward_integration_state,
                                                   0.0,
                                                   1.0,
                                                   m_dX);


        //  Save a track of the state
        m_forward_integration_state_X1 = forward_integration_state;
    }



    void updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                     const Eigen::VectorXd &t_Delta_dot_qe,
                                     const Eigen::VectorXd &t_Delta_ddot_qe)
    {
        m_Delta_qe      = t_Delta_qe;
        m_Delta_dot_qe  = t_Delta_dot_qe;
        m_Delta_ddot_qe = t_Delta_ddot_qe;
    }




    void forwardTangentKinematics()
    {
        forwardTangentKinematics(Eigen::Vector3d::Zero(),
                                 Eigen::Vector3d::Zero(),
                                 Eigen::Vector3d::Zero(),
                                 Eigen::Vector3d::Zero(),
                                 Eigen::Vector3d::Zero(),
                                 Eigen::Vector3d::Zero());
    }




    void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                  const Eigen::Vector3d &t_initial_Delta_position,
                                  const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                  const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                  const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                  const Eigen::Vector3d &t_initial_Delta_linear_acceleration)
    {
        ExplicitIntegrationODEs::TangentKinematicState tangent_kinematics_state;

        tangent_kinematics_state << m_forward_integration_state_X0,
                                    t_initial_Delta_orientation,
                                    t_initial_Delta_position,
                                    t_initial_Delta_angular_velocity,
                                    t_initial_Delta_linear_velocity,
                                    t_initial_Delta_angular_acceleration,
                                    t_initial_Delta_linear_acceleration;

        //  Save a track of the state
        m_forward_tangent_kinematics_state_X0 = tangent_kinematics_state;

    //  Forward Integration
        boost::numeric::odeint::integrate_adaptive(stepper,
                                                   [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                        {this->tangentKinematicsODEs(t_y, t_dyds, t_s);},
                                                   tangent_kinematics_state,
            0.0,
            1.0,
            m_dX);

        //  Save a track of the state
        m_forward_tangent_kinematics_state_X1 = tangent_kinematics_state;

    }




    ::LieAlgebra::Kinematics getKinematicsAtTip()const
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



    ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const
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




    void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1)
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

        Eigen::VectorXd Q1 = Eigen::VectorXd::Zero(m_generalised_coordinates_dimension);


        m_backward_integration_state_X1 <<
                m_forward_integration_state_X1,
                t_Lambda_X1,
                Q1;


        m_backward_integration_state_X0 = m_backward_integration_state_X1;



        //  Backward Integration
        boost::numeric::odeint::integrate_adaptive(stepper,
                                                   [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                    {this->backwardODEs(t_y, t_dyds, t_s);},
                                                    m_backward_integration_state_X0,
                                                    1.0,
                                                    0.0,
                                                    -m_dX);


    }






    void updateInternalActuation(const double &t_time)
    {
        m_distributed_actuation->updateActuation(t_time);

        m_Qad_X0.setZero();

        //  Backward Integration
        m_integrate_Qad();

    }


    Eigen::VectorXd getQad() const
    {
        return m_Qad_X0;
    }


    void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
    {

        /*  The state has the form
         *  | Q  |   w, x, y, z                       0-3
         *  | r  |   x, y, z                          4-6
         *  | η  |   Ω1, Ω2, Ω3, V1, V2, V3          7-12
         *  | η̇  |   Ω1, Ω2, Ω3, V1, V2, V3         13-18
         *  | Λ  |   C1, C2, C3, N1, N2, N3          19-24
         *  | ∆ζ |  ∆K1, ∆K2, ∆K3, ∆Γ1, ∆Γ2, ∆Γ3    25-30
         *  | ∆η |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   31-36
         *  | ∆η̇ |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   37-42
         *  | ∆Λ |   C1, C2, C3, N1, N2, N3          43-48
         *  | ∆Qa|                                   49-49+ne
         */

        Eigen::VectorXd Delta_Q1 = Eigen::VectorXd::Zero(m_generalised_coordinates_dimension);

        Eigen::VectorXd Lambda_X1 = m_backward_integration_state_X1.block<6, 1>(19,0);

        m_tangent_dynamics_state_X1 <<
                      m_forward_tangent_kinematics_state_X1.block<19,1>(0,0),
                      Lambda_X1,
                      m_forward_tangent_kinematics_state_X1.block<18,1>(19,0),
                      t_Delta_Lambda_X1,
                      Delta_Q1;


        m_tangent_dynamics_state_X0 = m_tangent_dynamics_state_X1;



    //  Backward Integration
        boost::numeric::odeint::integrate_adaptive(stepper,
                                                   [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                    {this->tangentDynamicsODEs(t_y, t_dyds, t_s);},
                                                          m_tangent_dynamics_state_X0,
                                                          1.0,
                                                          0.0,
                                                          -m_dX);

    }



    ::LieAlgebra::Vector6d getLambdaAtBase()const
    {
        return m_backward_integration_state_X0.block<6, 1>(19, 0);
    }



    LieAlgebra::Vector6d getDeltaLambdaAtBase()const
    {
        return m_tangent_dynamics_state_X0.block<6, 1>(43,0);
    }


    ::LieAlgebra::Vector6d getQaAtBase()const
    {
        return m_backward_integration_state_X0.block(25, 0, m_generalised_coordinates_dimension, 1);
    }



    LieAlgebra::Vector6d getDeltaQaAtBase()const
    {
        return m_tangent_dynamics_state_X0.block(49, 0, m_generalised_coordinates_dimension, 1);
    }


    void updateIntegrationDomain(const double &t_rod_lenght)
    {
        m_rod_length = t_rod_lenght;
    }






    Eigen::MatrixXd getRodPositions() const
    {

        /*  The forward state has the form
         *  | Q |   w, x, y, z                  0-3
         *  | r |   x, y, z                     4-6
         */
        Eigen::VectorXd y = m_forward_integration_state_X0.block<7, 1>(0, 0);


        Eigen::MatrixXd rod_shapes_stack { Eigen::MatrixXd::Zero(3, m_number_of_Chebyshev_points) };

    //  Forward Integration with observer
        boost::numeric::odeint::integrate_times(stepper,
                                                [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                 {this->forwardStaticODEs(t_y, t_dyds, t_s);},
        y,
        m_Chebyshev_points_reversed.begin(),
        m_Chebyshev_points_reversed.end(),
        m_dX,
        [&, col=0](const Eigen::VectorXd& t_y, const double)mutable{
            rod_shapes_stack.block<3, 1>(0, col++) = t_y.block<3,1>(4, 0);
        });


        return rod_shapes_stack;
    }




    FullODEStatesObservations getIDMStatesObservations()
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

        Eigen::VectorXd Q1 = Eigen::VectorXd::Zero(m_generalised_coordinates_dimension);



        m_backward_integration_state_X0 = m_backward_integration_state_X1;




        //  Backward Integration
        boost::numeric::odeint::integrate_times(stepper,
                                                [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                 {this->backwardODEs(t_y, t_dyds, t_s);},
            m_backward_integration_state_X0,
            m_Chebyshev_points.begin(),
            m_Chebyshev_points.end(),
            -m_dX,
        [this, col=m_number_of_Chebyshev_points](const Eigen::VectorXd& t_y, const double)mutable{
            m_IDM_states_observation.writeState(t_y, --col);
        });

        return m_IDM_states_observation;
    }





    FullODEStatesObservations getTIDMStatesObservations()
    {
        /*  The state has the form
         *  | Q  |   w, x, y, z                       0-3
         *  | r  |   x, y, z                          4-6
         *  | η  |   Ω1, Ω2, Ω3, V1, V2, V3          7-12
         *  | η̇  |   Ω1, Ω2, Ω3, V1, V2, V3         13-18
         *  | Λ  |   C1, C2, C3, N1, N2, N3          19-24
         *  | ∆ζ |  ∆K1, ∆K2, ∆K3, ∆Γ1, ∆Γ2, ∆Γ3    25-30
         *  | ∆η |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   31-36
         *  | ∆η̇ |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   37-42
         *  | ∆Λ |   C1, C2, C3, N1, N2, N3          43-48
         *  | ∆Qa|                                   49-49+ne
         */

        Eigen::VectorXd Delta_Q1 = Eigen::VectorXd::Zero(m_generalised_coordinates_dimension);

        Eigen::VectorXd Lambda_X1 = m_backward_integration_state_X1.block<6, 1>(19,0);

        m_tangent_dynamics_state_X0 = m_tangent_dynamics_state_X1;



        //  Backward Integration
        boost::numeric::odeint::integrate_times(stepper,
                                                [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                 {this->tangentDynamicsODEs(t_y, t_dyds, t_s);},
            m_tangent_dynamics_state_X0,
            m_Chebyshev_points.begin(),
            m_Chebyshev_points.end(),
            -m_dX,
        [this, col=m_number_of_Chebyshev_points](const Eigen::VectorXd& t_y, const double)mutable{
            m_TIDM_states_observation.writeState(t_y.block(25, 0, 24+m_generalised_coordinates_dimension, 1), --col);

        });

        return m_TIDM_states_observation;
    }



    Eigen::MatrixXd integratePhiTPhi()const
    {

        Eigen::VectorXd PhiTPhi_state = Eigen::VectorXd::Zero(m_generalised_coordinates_dimension*m_generalised_coordinates_dimension);



        boost::numeric::odeint::integrate_adaptive(stepper,
                                                   [this](const Eigen::VectorXd &t_y, Eigen::VectorXd &t_dyds, const double t_s)
                                                    {this->PhiTPhiODE(t_y, t_dyds, t_s);},
                                                    PhiTPhi_state,
                                                    0.0,
                                                    1.0,
                                                    m_dX);


        Eigen::MatrixXd PhiTPhi = Eigen::Map<Eigen::MatrixXd>(PhiTPhi_state.data(),
                                                              m_generalised_coordinates_dimension,
                                                              m_generalised_coordinates_dimension);


        return PhiTPhi;
    }


    unsigned int m_number_of_Chebyshev_points;

    double m_dX { 0.001 };


    Method stepper { Method() };




    std::vector<double> m_Chebyshev_points { ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points) };
    std::vector<double> m_Chebyshev_points_reversed { [this](){
            std::vector<double> Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points);

            std::vector<double> Chebyshev_points_reversed;


            for(int i=m_number_of_Chebyshev_points-1; i>=0; i--)
                Chebyshev_points_reversed.push_back( Chebyshev_points[i] );

            return Chebyshev_points_reversed;
        }()
    };

    //    Eigen::MatrixXd m_IDM_states_stack { Eigen::MatrixXd::Zero(m_backward_state_size, m_number_of_Chebyshev_points) };
    FullODEStatesObservations m_IDM_states_observation { FullODEStatesObservations(m_number_of_Chebyshev_points,
                                                                                   m_generalised_coordinates_dimension) };

    FullODEStatesObservations m_TIDM_states_observation { FullODEStatesObservations(m_number_of_Chebyshev_points,
                                                                                    m_generalised_coordinates_dimension,
                                                                                    {3, 1}) };


};




}   //  namespace CROSP::numerical_integrators::explicit_methods

#endif // EXPLICIT_METHOTDS_HPP
