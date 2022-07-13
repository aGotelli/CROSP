/**
 * \file integrators.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod IDM
 * \date 13-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */



#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <memory>

#include "OSNI/include/OSNI.hpp"


namespace CROSP {

struct QuaternionIntegrator : public OSNI::ODEA {
    QuaternionIntegrator() : OSNI::ODEA(4) {}

    QuaternionIntegrator(const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEA(4, t_number_of_Chebyshev_points) {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final
    {

    }

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;

};

}





#endif // INTEGRATORS_HPP
