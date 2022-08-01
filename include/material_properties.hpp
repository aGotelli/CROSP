/**
 * \file material_properties.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the structure defining the material properties of the rod
 * \date 1-08-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */



#ifndef MATERIAL_PROPERTIES_HPP
#define MATERIAL_PROPERTIES_HPP

#include <Eigen/Dense>
#include <math.h>


namespace CROSP {


struct MaterialProperties {

    MaterialProperties()=default;

    const double m_E { 210e9 };
    const double m_G {  80e9 };
    const double m_rho { 7800 };


    const double m_r { 0.01 };
    const double m_A { M_PI*m_r*m_r };
    const double m_L { 1.0 };

    const double m_Jx { M_PI*pow(m_r,4)/2 };
    const double m_Jy { M_PI*pow(m_r,4)/4 };
    const double m_Jz { M_PI*pow(m_r,4)/4 };




    const Eigen::Matrix<double, 6, 6> m_H { [&](){
            Eigen::Matrix<double, 6, 6> H;
            H.setZero();
            H.diagonal() << m_G*m_Jx, m_E*m_Jy, m_E*m_Jz, m_E*m_A, m_G*m_A, m_G*m_A;

            return H;}() };


    const  Eigen::Matrix<double, 6, 6>  m_M{ [&](){
            Eigen::Matrix<double, 6, 6> M = Eigen::Matrix<double, 6, 6>::Zero();

            M.diagonal() << m_rho*m_Jx, m_rho*m_Jy, m_rho*m_Jz, m_rho*m_A, m_rho*m_A, m_rho*m_A;
            return M;}() };

};








}

#endif // MATERIAL_PROPERTIES_HPP
