#include "CROSP/rod_properties/rod_properties.hpp"

#include <boost/numeric/odeint.hpp>

namespace CROSP::rod_properties {

RodDimensions::RodDimensions(CrossSection* t_cross_section, const double &t_L)
    : m_cross_section( t_cross_section ),
      m_L(t_L)
{}



RodDimensions::~RodDimensions()
{
    delete m_cross_section;
}






RodProperties::RodProperties(const std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation)
    : m_Kee( defineKee(t_strain_parameterisation->m_polynomial_representation) )
{}


RodProperties::RodProperties(const std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                             const RodDimensions &t_rod_dimensions,
                             const MaterialProperties &t_material_properties)
    : m_material_properties(t_material_properties),
      m_rod_dimensions(t_rod_dimensions),
      m_Kee( defineKee(t_strain_parameterisation->m_polynomial_representation) )
{}


void RodProperties::updateRodLength(const double &t_rod_lenght)
{
    const double scale = t_rod_lenght/m_rod_dimensions.m_L;
    m_rod_dimensions.m_L = t_rod_lenght;

    m_Kee *= scale;
    m_Dee *= scale;
}




double RodProperties::distributedDensity()const
{
    return m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Area();
}


Eigen::Vector3d RodProperties::distributedGravitationalForce()const
{
    return m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Area()*m_gravity;
}



Eigen::Matrix3d RodProperties::getMAngular()const
{
    return m_M.block<3, 3>(0, 0);
}

Eigen::Matrix3d RodProperties::getMLinear()const
{
    return m_M.block<3, 3>(3, 3);
}

Eigen::Matrix3d RodProperties::getHAngular()const
{
    return m_H.block<3, 3>(0, 0);
}

Eigen::Matrix3d RodProperties::getHLinear()const
{
    return m_H.block<3, 3>(3, 3);
}


Eigen::MatrixXd RodProperties::defineKee(const std::shared_ptr<const polynomial_representation::PolynomialRepresentation> t_polynomial_representation)const
{

    const unsigned int n = t_polynomial_representation->getCoordinatesDimension();

    const Eigen::MatrixXd Ha = t_polynomial_representation->m_B.transpose() * m_H * t_polynomial_representation->m_B;

    typedef Eigen::MatrixXd Kee_state_type;


    typedef boost::numeric::odeint::runge_kutta_dopri5< Kee_state_type, double,
                                                         Kee_state_type, double,
                                                         boost::numeric::odeint::vector_space_algebra> Ke_stepper;
    Eigen::MatrixXd Kee = Eigen::MatrixXd::Zero(n, n);

    boost::numeric::odeint::integrate_adaptive(Ke_stepper(), [&](const Kee_state_type &, Kee_state_type &t_dKeeds, const double t_X){
        const auto Phi = t_polynomial_representation->getPhi(t_X);

        t_dKeeds = Phi.transpose()*Ha*Phi;}, Kee, 0.0, m_rod_dimensions.m_L, 0.0005);

    return Kee;
}


}   //  namespace CROSP::rod_properties
