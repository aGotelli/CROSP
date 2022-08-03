#include "CROSP/material_properties/material_properties.hpp"

namespace CROSP::material_properties {


double MaterialProperties::distributedDensity()const
{
    return m_rho*m_A;
}

Eigen::Matrix3d MaterialProperties::getMAngular()const
{
    return m_M.block<3, 3>(0, 0);
}

Eigen::Matrix3d MaterialProperties::getMLinear()const
{
    return m_M.block<3, 3>(3, 3);
}

Eigen::Matrix3d MaterialProperties::getHAngular()const
{
    return m_H.block<3, 3>(0, 0);
}

Eigen::Matrix3d MaterialProperties::getHLinear()const
{
    return m_H.block<3, 3>(3, 3);
}


}   //  namespace CROSP::MaterialProperties
