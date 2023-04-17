/**
 * \file rod_properties.hpp
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


#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"


/// \brief CROSP::rod_properties namespace contains the definition of the rod properties
namespace CROSP::rod_properties {



/*!
 * \brief The MaterialProperties struct defines the material properties of the rod
 *
 * This structure contains the material properties of the rod. Namely the Young modulus E and the
 * shear modulus G, with the specific weight rho.
 */
struct MaterialProperties {

    MaterialProperties()=default;


    ///  \brief m_E Young modulus [Pa]
    double m_E { 210e9 };

    ///  \brief m_G Shear modulus [Pa]
    double m_G {  80e9 };

    ///  \brief m_rho Specific weight [kg/m^3]
    double m_rho { 7800 };

    /// \brief m_mu the internal dumping of the rod
    double m_mu { 1e-4 };

};


struct CrossSection{

    virtual ~CrossSection()=default;

    virtual double Area()const=0;

    virtual double Ixx()const=0;

    virtual double Iyy()const=0;

    virtual double Izz()const=0;

    virtual std::string printProperties() const=0;
};

struct CircularCrossSection : public CrossSection {

    CircularCrossSection(const double &t_radius=0.001)
        : m_radius(t_radius)
    {}

    virtual double Area()const final
    {
        return M_PI*m_radius*m_radius;
    }

    virtual double Ixx()const final
    {
        return M_PI*pow(m_radius,4)/2;
    }

    virtual double Iyy()const final
    {
        return M_PI*pow(m_radius,4)/4;
    }

    virtual double Izz()const final
    {
        return M_PI*pow(m_radius,4)/4;
    }

    virtual std::string printProperties() const final
    {
        std::stringstream properties;
        properties << "     r : " << m_radius << "\n";
        return properties.str();
    }

    double m_radius;
};


struct RectangularCrossSection : public CrossSection {

    RectangularCrossSection(const double &t_width,
                            const double &t_height)
        : m_width(t_width),
          m_height(t_height)
    {}

    virtual double Area()const final
    {
        return m_width*m_height;
    }

    virtual double Ixx()const final
    {
        return (m_width*m_height/12)*(m_width*m_width + m_height*m_height);
    }

    virtual double Iyy()const final
    {
        return m_width*pow(m_height, 3)/12;
    }

    virtual double Izz()const final
    {
        return pow(m_width, 3)*m_height/12;
    }

    virtual std::string printProperties() const final
    {
        std::stringstream properties;
        properties << "     width : " << m_width << ",";
        properties << "     height : " << m_height << "\n";
        return properties.str();
    }

    double m_width;
    double m_height;
};



/*!
 * \brief The RodDimensions struct contains the geometrical dimensions of the rod
 *
 * This object contains the geometrical foundamentals like radius and length.
 * It also computes the corresponding derived properties such as area and inertia modules.
 */
struct RodDimensions {

    RodDimensions()=default;

    ~RodDimensions();


    /*!
     * \brief RodDimensions construct the object given the properties
     * \param t_r   Radius of the section [m]
     * \param t_L   Length of the rod [m]
     */
    RodDimensions(CrossSection* t_cross_section, const double &t_L);


    CrossSection* m_cross_section {
        new CircularCrossSection()
    };

    /// \brief m_L  Length of the rod [m]
    double m_L { 1.0 };

};




/*!
 * \brief The RodProperties class contains the properties of the Cosserat rod
 */
class RodProperties {

public:


    RodProperties(const std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation);


    RodProperties(const std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                  const RodDimensions &t_rod_dimensions,
                  const MaterialProperties &t_material_properties);



    /// \brief m_material_properties instance of the rod material properties
    MaterialProperties m_material_properties;

    /// \brief m_rod_dimensions instance of the rod geometrical properties and dimensions
    RodDimensions m_rod_dimensions;


    /*!
     * \brief distributedDensity returns the distributed density of the rod
     * \return the distributed density of the rod
     */
    double distributedDensity()const;


    /*!
     * \brief distributedGravitationalForce returns the distributed force due to the weight of the rod
     * \return the distributed force due to the weight of the rod
     */
    Eigen::Vector3d distributedGravitationalForce()const;


    /*!
     * \brief getMAngular returns the angular subset of the inertia matrix
     * \return the angular subset of the inertia matrix
     */
    Eigen::Matrix3d getMAngular()const;


    /*!
     * \brief getMLinear returns the linear subset of the inertia matrix
     * \return the linear subset of the inertia matrix
     */
    Eigen::Matrix3d getMLinear()const;


    /*!
     * \brief getHAngular returns the angular subset of the Hookean matrix
     * \return the angular subset of the Hookean matrix
     */
    Eigen::Matrix3d getHAngular()const;


    /*!
     * \brief getHLinear returns the linear subset of the Hookean matrix
     * \return the linear subset of the Hookean matrix
     */
    Eigen::Matrix3d getHLinear()const;



    /// \brief m_H The Hookean matrix default initialised using the members m_material_properties and m_rod_dimensions
    const Eigen::Matrix<double, 6, 6> m_H { [&](){
            Eigen::Matrix<double, 6, 6> H;
            H.setZero();
            H.diagonal() << m_material_properties.m_G * m_rod_dimensions.m_cross_section->Ixx(),
                            m_material_properties.m_E * m_rod_dimensions.m_cross_section->Iyy(),
                            m_material_properties.m_E * m_rod_dimensions.m_cross_section->Izz(),
                            m_material_properties.m_E * m_rod_dimensions.m_cross_section->Area(),
                            m_material_properties.m_G * m_rod_dimensions.m_cross_section->Area(),
                            m_material_properties.m_G * m_rod_dimensions.m_cross_section->Area();

            return H;}() };

    /// \brief m_M The inertia matrix default initialised using the members m_material_properties and m_rod_dimensions
    const  Eigen::Matrix<double, 6, 6>  m_M{ [&](){
            Eigen::Matrix<double, 6, 6> M = Eigen::Matrix<double, 6, 6>::Zero();

            M.diagonal() << m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Ixx(),
                            m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Iyy(),
                            m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Izz(),
                            m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Area(),
                            m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Area(),
                            m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Area();
            return M;}() };


    /// \brief m_Kee The generalised elasticity matrix
    const Eigen::MatrixXd m_Kee;

    /// \brief m_Dee The matrix of the internal dumping
    const Eigen::MatrixXd m_Dee { m_material_properties.m_mu*m_Kee };


private:

    /*!
     * \brief defineKee defines the elasticity matrix Kee
     * \param t_ne the number of modes per admitted deformation
     * \param t_na the number of deformations degrees of freedom
     * \param t_B the map matrix to map the allowed strains in the space of the full strain
     * \return
     */
    Eigen::MatrixXd defineKee(const std::shared_ptr<const polynomial_representation::PolynomialRepresentation> t_polynomial_representation)const;


    const double gamma = 9.81;
public : Eigen::Vector3d m_gravity { Eigen::Vector3d(0, 0, -gamma) };

};







}   //  namespace CROSP::rod_properties

#endif // ROD_PROPERTIES_HPP
