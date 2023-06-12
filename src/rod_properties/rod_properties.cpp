#include "CROSP/rod_properties/rod_properties.hpp"

#include <boost/numeric/odeint.hpp>

namespace CROSP::rod_properties {






RodProperties::RodProperties(polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                             RodDimensions t_rod_dimensions,
                             MaterialProperties t_material_properties)
    : m_material_properties(t_material_properties),
      m_rod_dimensions( RodDimensions(std::move(t_rod_dimensions.m_cross_section),
                                      t_rod_dimensions.m_L)),
      m_Kee( defineKee( t_polynomial_representation ) )
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




void RodProperties::updateRodProperties(const double &t_EI,
                                        const polynomial_representation::PolynomialRepresentation &t_polynomial_representation)
{
    m_H(2, 2) = t_EI;

    m_Kee = defineKee(t_polynomial_representation);

    m_Dee = m_material_properties.m_mu*m_Kee;
}


::LieAlgebra::Matrix6d RodProperties::computeHookTensorMatrix()const
{
    ::LieAlgebra::Matrix6d H = ::LieAlgebra::Matrix6d::Zero();
    H.diagonal() << m_material_properties.m_G * m_rod_dimensions.m_cross_section->Ixx(),
                    m_material_properties.m_E * m_rod_dimensions.m_cross_section->Iyy(),
                    m_material_properties.m_E * m_rod_dimensions.m_cross_section->Izz(),
                    m_material_properties.m_E * m_rod_dimensions.m_cross_section->Area(),
                    m_material_properties.m_G * m_rod_dimensions.m_cross_section->Area(),
                    m_material_properties.m_G * m_rod_dimensions.m_cross_section->Area();

    return H;
}

::LieAlgebra::Matrix6d RodProperties::computeCrossSectionalInertiaMatrix()const
{
    ::LieAlgebra::Matrix6d M = ::LieAlgebra::Matrix6d::Zero();
    M.diagonal() << m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Ixx(),
                    m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Iyy(),
                    m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Izz(),
                    m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Area(),
                    m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Area(),
                    m_material_properties.m_rho * m_rod_dimensions.m_cross_section->Area();
    return M;
}

Eigen::MatrixXd RodProperties::defineKee(const polynomial_representation::PolynomialRepresentation &t_polynomial_representation)const
{

    const unsigned int n = t_polynomial_representation.getCoordinatesDimension();

    const Eigen::MatrixXd Ha = t_polynomial_representation.m_B.transpose() * m_H * t_polynomial_representation.m_B;


    typedef boost::numeric::odeint::runge_kutta_dopri5< Eigen::MatrixXd, double,
                                                         Eigen::MatrixXd, double,
                                                         boost::numeric::odeint::vector_space_algebra> Ke_stepper;
    Eigen::MatrixXd Kee = Eigen::MatrixXd::Zero(n, n);

    const double X0 = 0.0;
    const double X1 = 1.0;
    const double dX = 0.0005;

    boost::numeric::odeint::integrate_adaptive(Ke_stepper(), [&](const Eigen::MatrixXd &, Eigen::MatrixXd &t_dKeeds, const double t_X){
        const auto Phi = t_polynomial_representation.getPhi( t_X );

        t_dKeeds = Phi.transpose()*Ha*Phi;
    }, Kee, X0, X1, dX);

    Kee *= m_rod_dimensions.m_L;

    return Kee;
}


}   //  namespace CROSP::rod_properties
