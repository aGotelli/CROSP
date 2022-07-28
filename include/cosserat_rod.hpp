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


#include <boost/math/special_functions/legendre.hpp>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "integrators.hpp"
#include "strain_parameterisation.hpp"

namespace CROSP {


struct MaterialProperties {

    const double m_E { 210000000000 };
    const double m_G {  80000000000 };
    const double m_rho { 7800 };


    const double m_r { 0.001 };
    const double m_A { M_PI*m_r*m_r };
    const double m_L { 1.0 };

    const double m_Jx { M_PI*pow(m_r,4)/2 };
    const double m_Jy { M_PI*pow(m_r,4)/4 };
    const double m_Jz { M_PI*pow(m_r,4)/4 };




    const Eigen::Matrix<double, 6, 6> m_H { [this](){
            Eigen::Matrix<double, 6, 6> H;
            H.setZero();
            H.diagonal() << m_G*m_Jx, m_E*m_Jy, m_E*m_Jz, m_E*m_A, m_G*m_A, m_G*m_A;

            return H;}() };


    const  Eigen::Matrix<double, 6, 6>  m_M{ Eigen::Matrix<double, 6, 6>::Identity()/*[this](){
            Eigen::Matrix<double, 6, 6> M = Eigen::Matrix<double, 6, 6>::Zero();

            M.diagonal() << m_rho*m_Jx, m_rho*m_Jy, m_rho*m_Jz, m_rho*m_A, m_rho*m_A, m_rho*m_A;
            return M;}()*/ };

};


class CosseratRod
{
public:
    CosseratRod()=default;

    CosseratRod(const unsigned int t_number_of_chebyshev_points);

    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                const Eigen::VectorXd &t_dot_qe,
                                const Eigen::VectorXd &t_ddot_qe);

    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion=Eigen::Vector4d(1, 0, 0, 0),
                           const Eigen::Vector3d &t_initial_position=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_angular_velocity=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_linear_velocity=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_angular_acceleration=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_linear_acceleration=Eigen::Vector3d::Zero());

    ::LieAlgebra::Kinematics getKinematicsAtTip()const;



    void backwardDynamics(const Eigen::Vector3d &t_force_at_tip,
                          const Eigen::Vector3d &t_couple_at_tip);



#ifndef DEVELOPER
private:
#endif


    const unsigned int m_ne { 3 };
    const unsigned int m_na { 3 };


    const unsigned int m_number_of_chebyshev_points { 17 };

    const std::vector<double> m_Chebyshev_points { ::Chebyshev::ComputeChebyshevPoints(m_number_of_chebyshev_points) };

    const std::function<double(const unsigned int, const double&)> m_polynomial_base { [](const unsigned int t_point, const double& t_x) {return boost::math::legendre_p(t_point, t_x);} };

    MaterialProperties m_material_properties;

    StrainParameterisation m_strain_parameterisation { StrainParameterisation(generatePhiStack(m_ne, m_na, m_Chebyshev_points, m_polynomial_base)) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Lambda_stack{ std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_chebyshev_points) };



    std::shared_ptr<OSNI::ODESolverInterface> m_quaternion_integrator { std::make_shared<QuaternionIntegrator>(m_K_stack,
                                                                                                               m_number_of_chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_position_integrator { std::make_shared<PositionIntegrator>(m_quaternion_integrator,
                                                                                                           m_Lambda_stack,
                                                                                                           m_number_of_chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_angular_velocity_integrator { std::make_shared<AngularVelocityIntegrator>(m_K_stack,
                                                                                                                          m_dot_K_stack,
                                                                                                                          m_number_of_chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_linear_velocity_integrator { std::make_shared<LinearVelocityIntegrator>(m_K_stack,
                                                                                                                        m_Lambda_stack,
                                                                                                                        m_dot_Lambda_stack,
                                                                                                                        m_angular_velocity_integrator,
                                                                                                                        m_number_of_chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_angular_acceleration_integrator { std::make_shared<AngularAccelerationIntegrator>(m_K_stack,
                                                                                                                                  m_dot_K_stack,
                                                                                                                                  m_ddot_K_stack,
                                                                                                                                  m_angular_velocity_integrator,
                                                                                                                                  m_number_of_chebyshev_points ) };

    std::shared_ptr<OSNI::ODESolverInterface> m_linear_acceleration_integrator { std::make_shared<LinearAccelerationIntegrator>(m_K_stack,
                                                                                                                                m_dot_K_stack,
                                                                                                                                m_Lambda_stack,
                                                                                                                                m_dot_Lambda_stack,
                                                                                                                                m_ddot_Lambda_stack,
                                                                                                                                m_angular_velocity_integrator,
                                                                                                                                m_linear_velocity_integrator,
                                                                                                                                m_angular_acceleration_integrator,
                                                                                                                                m_number_of_chebyshev_points ) };


    std::shared_ptr<OSNI::ODESolverInterface> m_internal_forces_integrator { std::make_shared<InternalForcesIntegrator>(m_material_properties.m_M.block<3,3>(3, 3),
                                                                                                                        m_K_stack, m_angular_velocity_integrator,
                                                                                                                        m_linear_velocity_integrator,
                                                                                                                        m_linear_acceleration_integrator,
                                                                                                                        m_quaternion_integrator,
                                                                                                                        m_position_integrator,
                                                                                                                        m_number_of_chebyshev_points )};

    std::shared_ptr<OSNI::ODESolverInterface> m_internal_couples_integrator { std::make_shared<InternalCouplesIntegrator>(m_material_properties.m_M.block<3,3>(0, 0),
                                                                                                                          m_material_properties.m_M.block<3,3>(3, 3),
                                                                                                                          m_K_stack,
                                                                                                                          m_Lambda_stack,
                                                                                                                          m_angular_velocity_integrator,
                                                                                                                          m_linear_velocity_integrator,
                                                                                                                          m_angular_acceleration_integrator,
                                                                                                                          m_quaternion_integrator,
                                                                                                                          m_position_integrator,
                                                                                                                          m_internal_forces_integrator,
                                                                                                                          m_number_of_chebyshev_points )};

};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
