#include "CROSP/rod_properties/rod_properties.hpp"

#include <boost/numeric/odeint.hpp>

namespace CROSP::rod_properties {

RodProperties::RodProperties(const double t_mu)
    : m_mu( t_mu )
{}



RodProperties::RodProperties(const double t_mu,
                             Eigen::Vector3d t_gravity)
    : m_mu( t_mu ),
      m_gravity(t_gravity)
{}


RodProperties::RodProperties(const double t_mu,
                             Eigen::Vector3d t_gravity,
                             const MaterialProperties &t_material_properties)
    : m_material_properties(t_material_properties),
      m_mu( t_mu ),
      m_gravity(t_gravity)
{}


RodProperties::RodProperties(const double t_mu,
                             const polynomial_representation::PolynomialRepresentation t_polynomial_representation)
    : m_Kee( defineKee(t_polynomial_representation) ),
      m_mu( t_mu )
{}


RodProperties::RodProperties(const polynomial_representation::PolynomialRepresentation t_polynomial_representation)
    : m_Kee( defineKee(t_polynomial_representation) )
{}


RodProperties::RodProperties(const MaterialProperties &t_material_properties)
    : m_material_properties(t_material_properties)
{}


RodProperties::RodProperties(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                             const MaterialProperties &t_material_properties)
    : m_material_properties(t_material_properties),
      m_Kee( defineKee(t_polynomial_representation) )
{}


RodProperties::RodProperties(const RodDimensions &t_rod_dimensions)
    : m_rod_dimensions(t_rod_dimensions)
{}

RodProperties::RodProperties(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                             const RodDimensions &t_rod_dimensions)
    : m_rod_dimensions(t_rod_dimensions),
      m_Kee( defineKee(t_polynomial_representation) )
{}

RodProperties::RodProperties(const MaterialProperties &t_material_properties,
                             const RodDimensions &t_rod_dimensions)
    : m_material_properties(t_material_properties),
      m_rod_dimensions(t_rod_dimensions)
{}


RodProperties::RodProperties(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                             const MaterialProperties &t_material_properties,
                             const RodDimensions &t_rod_dimensions)
    : m_material_properties(t_material_properties),
      m_rod_dimensions(t_rod_dimensions),
      m_Kee( defineKee(t_polynomial_representation) )
{}




double RodProperties::distributedDensity()const
{
    return m_material_properties.m_rho * m_rod_dimensions.m_A;
}


Eigen::Vector3d RodProperties::distributedGravitationalForce()const
{
    return m_material_properties.m_rho * m_rod_dimensions.m_A*m_gravity;
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


Eigen::MatrixXd RodProperties::defineKee(const polynomial_representation::PolynomialRepresentation &t_polynomial_representation)const
{

    const unsigned int n = t_polynomial_representation.getCoordinatesDimension();

    const Eigen::MatrixXd Ha = t_polynomial_representation.m_B.transpose() * m_H * t_polynomial_representation.m_B;

    typedef Eigen::MatrixXd Kee_state_type;


    typedef boost::numeric::odeint::runge_kutta_dopri5< Kee_state_type, double,
                                                         Kee_state_type, double,
                                                         boost::numeric::odeint::vector_space_algebra> Ke_stepper;
    Eigen::MatrixXd Kee = Eigen::MatrixXd::Zero(n, n);

    boost::numeric::odeint::integrate_adaptive(Ke_stepper(), [&](const Kee_state_type &, Kee_state_type &t_dKeeds, const double t_X){
        const auto Phi = t_polynomial_representation.getPhi(t_X);

        t_dKeeds = Phi.transpose()*Ha*Phi;}, Kee, 0.0, 1.0, 0.0005);

    return Kee;
}


}   //  namespace CROSP::rod_properties
