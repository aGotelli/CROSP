#include "CROSP/rod_properties/rod_properties.hpp"

namespace CROSP::rod_properties {


RodProperties::RodProperties(const MaterialProperties &t_material_properties) :
    m_material_properties(t_material_properties)
{}

RodProperties::RodProperties(const RodDimensions &t_rod_dimensions) :
    m_rod_dimensions(t_rod_dimensions)
{}

RodProperties::RodProperties(const MaterialProperties &t_material_properties,
                             const RodDimensions &t_rod_dimensions) :
    m_material_properties(t_material_properties),
    m_rod_dimensions(t_rod_dimensions)
{}



double RodProperties::distributedDensity()const
{
    return m_material_properties.m_rho * m_rod_dimensions.m_A;
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


}   //  namespace CROSP::rod_properties
