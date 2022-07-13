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



#ifndef COSSERAT_ROD_HPP
#define COSSERAT_ROD_HPP

#include <Eigen/Dense>

#include "integrators.hpp"

namespace CROSP {


class CosseratRod
{
public:
    CosseratRod();

    void updateParameterisation(const Eigen::VectorXd t_q,
                                const Eigen::VectorXd t_dot_q,
                                const Eigen::VectorXd t_ddot_q);

private:


    const unsigned int m_number_of_chebyshev_points { 16 };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack;
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Lambda_stack;

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack;
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack;
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack;


    std::shared_ptr<OSNI::ODESolverInterface> m_quaternion_integrator;

    std::shared_ptr<OSNI::ODESolverInterface> m_position_integrator;

};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
