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



#ifndef ROD_PROPERTIES_HPP
#define ROD_PROPERTIES_HPP

#include <Eigen/Dense>
#include <math.h>


#include "CROSP/base_maps/base_maps.hpp"


namespace CROSP::rod_properties {



/*!
 * \brief The MaterialProperties struct defines the material properties of the rod
 *
 * This structure contains the material properties of the rod. Namely the Young modulus E and the
 * shear modulus G, with the specific weight rho.
 */
struct MaterialProperties {

    MaterialProperties()=default;


    /*!
     * \brief MaterialProperties construct the object with the given properties
     * \param t_E   Young modulus [Pa]
     * \param t_G   Shear modulus [Pa]
     * \param t_rho Specific weight [kg/m^3]
     */
    MaterialProperties(const double &t_E, const double &t_G, const double &t_rho) : m_E(t_E),
                                                                                    m_G(t_G),
                                                                                    m_rho(t_rho)
    {}

    //  Young modulus [Pa]
    const double m_E { 210e9 };

    //  Shear modulus [Pa]
    const double m_G {  80e9 };

    //  Specific weight [kg/m^3]
    const double m_rho { 7800 };

};



/*!
 * \brief The RodDimensions struct contains the geometrical dimensions of the rod
 *
 * This object contains the geometrical foundamentals like radius and length.
 * It also computes the corresponding derived properties such as area and inertia modules.
 */
struct RodDimensions {

    RodDimensions()=default;

    /*!
     * \brief RodDimensions construct the object given the properties
     * \param t_r   Radius of the section [m]
     * \param t_L   Length of the rod [m]
     */
    RodDimensions(const double &t_r, const double &t_L) : m_r(t_r),
                                                          m_L(t_L)
    {}


    //  Radius of the section [m]
    const double m_r { 0.01 };

    //  Area of the section [m^2]
    const double m_A { M_PI*m_r*m_r };

    //  Length of the rod [m]
    const double m_L { 1.0 };

    //  Moment of inertia [m^4]
    const double m_Jx { M_PI*pow(m_r,4)/2 };
    const double m_Jy { M_PI*pow(m_r,4)/4 };
    const double m_Jz { M_PI*pow(m_r,4)/4 };
};





class RodProperties {

public:

    RodProperties()=default;

    RodProperties(const unsigned int t_ne,
                  const unsigned int t_na,
                  const Eigen::MatrixXd &t_B);

    RodProperties(const MaterialProperties &t_material_properties);

    RodProperties(const RodDimensions &t_rod_dimensions);

    RodProperties(const MaterialProperties &t_material_properties,
                  const RodDimensions &t_rod_dimensions);

    RodProperties(const unsigned int t_ne,
                  const unsigned int t_na,
                  const Eigen::MatrixXd &t_B,
                  const MaterialProperties &t_material_properties,
                  const RodDimensions &t_rod_dimensions);



    const MaterialProperties m_material_properties { MaterialProperties() };

    const RodDimensions m_rod_dimensions { RodDimensions() };


    double distributedDensity()const;

    Eigen::Matrix3d getMAngular()const;

    Eigen::Matrix3d getMLinear()const;

    Eigen::Matrix3d getHAngular()const;

    Eigen::Matrix3d getHLinear()const;




    const Eigen::Matrix<double, 6, 6> m_H { [&](){
            Eigen::Matrix<double, 6, 6> H;
            H.setZero();
            H.diagonal() << m_material_properties.m_G * m_rod_dimensions.m_Jx,
                            m_material_properties.m_E * m_rod_dimensions.m_Jy,
                            m_material_properties.m_E * m_rod_dimensions.m_Jz,
                            m_material_properties.m_E * m_rod_dimensions.m_A,
                            m_material_properties.m_G * m_rod_dimensions.m_A,
                            m_material_properties.m_G * m_rod_dimensions.m_A;

            return H;}() };


    const  Eigen::Matrix<double, 6, 6>  m_M{ [&](){
            Eigen::Matrix<double, 6, 6> M = Eigen::Matrix<double, 6, 6>::Zero();

            M.diagonal() << m_material_properties.m_rho * m_rod_dimensions.m_Jx,
                            m_material_properties.m_rho * m_rod_dimensions.m_Jy,
                            m_material_properties.m_rho * m_rod_dimensions.m_Jz,
                            m_material_properties.m_rho * m_rod_dimensions.m_A,
                            m_material_properties.m_rho * m_rod_dimensions.m_A,
                            m_material_properties.m_rho * m_rod_dimensions.m_A;
            return M;}() };


    const std::array<bool, 6> m_admitted_deformations { true, true, true, false, false, false };

    const unsigned int m_na { [&]()->unsigned int{ return std::count(m_admitted_deformations.begin(), m_admitted_deformations.end(), true);}() };

    const unsigned int m_ne { 3 };

    const Eigen::MatrixXd m_Kee { defineKee( m_ne, m_na, base_maps::getB(m_admitted_deformations) ) };


private:

    Eigen::MatrixXd defineKee(const unsigned int t_ne,
                              const unsigned int t_na,
                              const Eigen::MatrixXd &t_B)const;

};







}   //  namespace CROSP::rod_properties

#endif // ROD_PROPERTIES_HPP
