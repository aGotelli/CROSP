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

#include "idm_integrators.hpp"
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

    MaterialProperties m_material_properties;

    std::shared_ptr<StrainParameterisation> m_strain_parameterisation { std::make_shared<StrainParameterisation>() };

    std::shared_ptr<IDMIntegrators> m_idm_integrators { std::make_shared<IDMIntegrators>(m_strain_parameterisation, m_material_properties.m_M.block<3,3>(0,0), m_material_properties.m_M.block<3,3>(3,3))};

};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
