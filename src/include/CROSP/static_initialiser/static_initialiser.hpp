/**
 * \file static_initialiser.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains class for the static initialiser of a Cosserat rod
 * \date 15-09-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#ifndef STATIC_INITIALISER_HPP
#define STATIC_INITIALISER_HPP

#include "CROSP/CROSP/cosserat_rod.hpp"
#include <ceres/ceres.h>


namespace CROSP::static_initialiser {


class Initializer {
public:
    Initializer(::CROSP::CosseratRod &t_rod,
                ::LieAlgebra::Vector6d &t_F1,
                Eigen::VectorXd &t_initial_guess,
                bool init_google_logging=true);

    bool operator()(const double* const* t_guess, double* t_residual) const;

private:

    typedef ceres::DynamicNumericDiffCostFunction<Initializer,
                                                  ceres::NumericDiffMethodType::FORWARD> CostFunction;

    CostFunction* m_cost_function;

    std::shared_ptr<::CROSP::CosseratRod> m_rod_ptr { nullptr };

    const Eigen::VectorXd m_F1 { Eigen::VectorXd::Zero(6) };

    const unsigned int m_n { m_rod_ptr->getCoordinatesDimension() };

    //  Define the problem
    ceres::Problem m_problem;

    //  Set solver options
    ceres::Solver::Options m_options;

    //  Solver verbosity
    ceres::Solver::Summary m_summary;



    const Eigen::Quaterniond m_Q0;
    const Eigen::Vector3d m_r0;
};

}   //  namespace CROSP::static_initialiser

#endif // STATIC_INITIALISER_HPP
