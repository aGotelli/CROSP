#ifndef RUNGEKUTTAINTEGRATOR_HPP
#define RUNGEKUTTAINTEGRATOR_HPP



#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"



namespace CROSP::numerical_integrators::runge_kutta {

class RungeKuttaIntegrator : public CosseratRodIntegrators
{
public:
    RungeKuttaIntegrator(const strain_parameterisation::StrainParameterisation t_strain_parameterisation,
                         const strain_parameterisation::StrainParameterisation t_strain_parameterisation_Delta,
                         const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                         const rod_properties::RodProperties t_rod_properties,
                         const unsigned int t_number_of_Chebyshev_points);


    virtual std::string printIntegratorProperties()const final
    {
        std::stringstream integrator_properties;
        integrator_properties << "Runge Kutta integrator\n" "Adaptive step with initial value of : " << m_dX << " m";

        return integrator_properties.str();

    }


    virtual void updateParameterisation(const Eigen::VectorXd &t_qe,
                                        const Eigen::VectorXd &t_dot_qe,
                                        const Eigen::VectorXd &t_ddot_qe) final
    {
        m_qe      = t_qe;
        m_dot_qe  = t_dot_qe;
        m_ddot_qe = t_ddot_qe;
    }


    virtual void forwardKinematics() final;


    virtual void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                   const Eigen::Vector3d &t_initial_position,
                                   const Eigen::Vector3d &t_initial_angular_velocity,
                                   const Eigen::Vector3d &t_initial_linear_velocity,
                                   const Eigen::Vector3d &t_initial_angular_acceleration,
                                   const Eigen::Vector3d &t_initial_linear_acceleration) final;


    virtual void updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                             const Eigen::VectorXd &t_Delta_dot_qe,
                                             const Eigen::VectorXd &t_Delta_ddot_qe) final
    {
        m_Delta_qe      = t_Delta_qe;
        m_Delta_dot_qe  = t_Delta_dot_qe;
        m_Delta_ddot_qe = t_Delta_ddot_qe;
    }

    virtual void forwardTangentKinematics() final
    {
//        //  Integrate Delta zeta
//        m_tidm_integrators->m_Delta_rotation->solveSystem();
//        m_tidm_integrators->m_Delta_position->solveSystem();

//        //  Integrate Delta eta
//        m_tidm_integrators->m_Delta_angular_velocity->solveSystem();
//        m_tidm_integrators->m_Delta_linear_velocity->solveSystem();

//        //  Integrate Delta dot eta
//        m_tidm_integrators->m_Delta_angular_acceleration->solveSystem();
//        m_tidm_integrators->m_Delta_linear_acceleration->solveSystem();
    }


    virtual void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                          const Eigen::Vector3d &t_initial_Delta_position,
                                          const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                          const Eigen::Vector3d &t_initial_Delta_linear_acceleration) final
    {
//        //  Integrate Delta zeta
//        m_tidm_integrators->m_Delta_rotation->integrate( t_initial_Delta_orientation );
//        m_tidm_integrators->m_Delta_position->integrate( t_initial_Delta_position );

//        //  Integrate Delta eta
//        m_tidm_integrators->m_Delta_angular_velocity->integrate( t_initial_Delta_angular_velocity );
//        m_tidm_integrators->m_Delta_linear_velocity->integrate( t_initial_Delta_linear_velocity );

//        //  Integrate Delta dot eta
//        m_tidm_integrators->m_Delta_angular_acceleration->integrate( t_initial_Delta_angular_acceleration );
//        m_tidm_integrators->m_Delta_linear_acceleration->integrate( t_initial_Delta_linear_acceleration );

    }






    virtual ::LieAlgebra::Kinematics getKinematicsAtTip()const final
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


    virtual ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const final
    {
//        ::LieAlgebra::TangentKinematics rod_tip_tangent_kinematics;

//        rod_tip_tangent_kinematics.m_Delta_zeta <<
//                m_tidm_integrators->m_Delta_rotation->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
//                m_tidm_integrators->m_Delta_position->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

//        rod_tip_tangent_kinematics.m_Delta_twist <<
//                m_tidm_integrators->m_Delta_angular_velocity->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
//                m_tidm_integrators->m_Delta_linear_velocity->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

//        rod_tip_tangent_kinematics.m_Delta_acceleration <<
//                m_tidm_integrators->m_Delta_angular_acceleration->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END),
//                m_tidm_integrators->m_Delta_linear_acceleration->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

//        return rod_tip_tangent_kinematics;
    }



    virtual void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1) final;



    virtual void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1) final
    {

//        Eigen::Vector3d Delta_couple_at_tip = t_Delta_Lambda_X1.block<3, 1>(0, 0);
//        Eigen::Vector3d Delta_force_at_tip  = t_Delta_Lambda_X1.block<3, 1>(3, 0);


//        m_tidm_integrators->m_Delta_internal_forces->integrate(Delta_force_at_tip);
//        m_tidm_integrators->m_Delta_internal_couples->integrate(Delta_couple_at_tip);

//        m_tidm_integrators->m_Delta_generalised_forces->solveSystem();

    }




    virtual ::LieAlgebra::Vector6d getLambdaAtBase()const final
    {
        return m_backward_integration_state_X0.block<6, 1>(19, 0);
    }



    virtual LieAlgebra::Vector6d getDeltaLambdaAtBase()const final
    {
//        ::LieAlgebra::Vector6d Delta_Lambda;
//        Delta_Lambda <<   m_tidm_integrators->m_Delta_internal_couples->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN),
//                          m_tidm_integrators->m_Delta_internal_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

//        return Delta_Lambda;
    }


    virtual ::LieAlgebra::Vector6d getQaAtBase()const final
    {
        return m_backward_integration_state_X0.block(25, 0, m_ne, 1);
    }



    virtual LieAlgebra::Vector6d getDeltaQaAtBase()const final
    {
//        return m_tidm_integrators->m_Delta_generalised_forces->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);;
    }


    virtual void updateIntegrationDomain(const double &t_rod_lenght) final
    {
        m_rod_length = t_rod_lenght;
    }


    virtual Eigen::MatrixXd getRodPositions()const;





private:

    strain_parameterisation::StrainParameterisation m_strain_parameterisation;
    strain_parameterisation::StrainParameterisation m_strain_parameterisation_Delta;
    polynomial_representation::PolynomialRepresentation m_polynomial_representation;
    rod_properties::RodProperties m_rod_properties;

    double m_dX { 0.005 };

    unsigned int m_ne { m_polynomial_representation.getCoordinatesDimension() };


    double m_rod_length { 1.0 };


    unsigned int m_number_of_Chebyshev_points { 31 };
    std::vector<double> m_Chebyshev_points_reversed { [this](){
            std::vector<double> Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points);

            std::vector<double> Chebyshev_points_reversed;


            for(int i=m_number_of_Chebyshev_points-1; i>=0; i--)
                Chebyshev_points_reversed.push_back( Chebyshev_points[i] );

            return Chebyshev_points_reversed;
        }()
    };





    Eigen::VectorXd m_forward_integration_state_X0 { Eigen::VectorXd::Zero(19) };
    Eigen::VectorXd m_forward_integration_state_X1 { Eigen::VectorXd::Zero(19) };

    Eigen::VectorXd m_backward_integration_state_X0 { Eigen::VectorXd::Zero(25+m_ne) };


    Eigen::VectorXd m_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };
    Eigen::VectorXd m_dot_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };
    Eigen::VectorXd m_ddot_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };


    Eigen::VectorXd m_Delta_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };
    Eigen::VectorXd m_Delta_dot_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };
    Eigen::VectorXd m_Delta_ddot_qe { Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension()) };




    void forwardODEs(const Eigen::VectorXd &t_y,
                     Eigen::VectorXd &t_dyds,
                     const double t_s) const;


    void backwardODEs(const Eigen::VectorXd &t_y,
                      Eigen::VectorXd &t_dyds,
                      const double t_s)const;



};

}   //  namespace CROSP::numerical_integrators::runge_kutta

#endif // RUNGEKUTTAINTEGRATOR_HPP
