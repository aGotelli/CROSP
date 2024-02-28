#include "CROSP/rod_properties/rod_properties.hpp"

#include <boost/numeric/odeint.hpp>

namespace CROSP::rod_properties {






RodProperties::RodProperties(RodDimensions t_rod_dimensions,
                             MaterialProperties t_material_properties)
    : m_material_properties(t_material_properties),
      m_rod_dimensions( RodDimensions(std::move(t_rod_dimensions.m_cross_section),
                                      t_rod_dimensions.m_L))
{}







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



}   //  namespace CROSP::rod_properties
