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

    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                const Eigen::VectorXd &t_dot_qe,
                                const Eigen::VectorXd &t_ddot_qe);

    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion=Eigen::Vector4d(1, 0, 0, 0),
                           const Eigen::Vector3d &t_initial_position=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_angular_velocity=Eigen::Vector3d::Zero());

private:



    const unsigned int m_ne { 3 };
    const unsigned int m_na { 3 };


    const unsigned int m_number_of_chebyshev_points { 16 };

    const std::vector<double> m_Chebyshev_points;

    const std::function<double(const unsigned int, const double&)> m_polynomial_base;

    Eigen::MatrixXd getPhi(const double& t_X,
                            const double& t_begin=0,
                            const double& t_end=1) const;




    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack;
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Lambda_stack;

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack;
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack;
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack;



    std::shared_ptr<OSNI::ODESolverInterface> m_quaternion_integrator;

    std::shared_ptr<OSNI::ODESolverInterface> m_position_integrator;

    std::shared_ptr<OSNI::ODESolverInterface> m_angular_velocity_integrator;

};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
