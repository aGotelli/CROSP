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

#include "cosserat_rod.hpp"

namespace CROSP {

CosseratRod::CosseratRod() :
    m_quaternion_integrator(std::make_shared<QuaternionIntegrator>(m_K_stack)),
    m_position_integrator(std::make_shared<PositionIntegrator>(m_quaternion_integrator, m_Lambda_stack))
{}



}   //  namespace CROSP
