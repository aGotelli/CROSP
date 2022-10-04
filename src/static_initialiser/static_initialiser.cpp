/**
 * \file static_initialiser.cpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains class for the static initialiser of a Cosserat rod
 * \date 15-09-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */
#include "CROSP/static_initialiser/static_initialiser.hpp"


namespace CROSP::static_initialiser {


Initializer::Initializer(::CROSP::CosseratRod &t_rod,
                         ::LieAlgebra::Vector6d &t_F1,
                         Eigen::VectorXd &t_initial_guess,
                         bool init_google_logging)
        : m_rod_ptr(std::make_shared<::CROSP::CosseratRod>(t_rod)),
          m_F1(t_F1)
{
    if(init_google_logging)
        google::InitGoogleLogging(std::string().c_str());

    CostFunction* cost_function = new CostFunction( this, ceres::DO_NOT_TAKE_OWNERSHIP );

    m_options.max_num_iterations = 100;
    m_options.linear_solver_type = ceres::DENSE_QR;
    m_options.minimizer_progress_to_stdout = true;

    //  Add the cost function to the problem
    m_problem.AddResidualBlock(cost_function, nullptr, t_initial_guess.data());

    //  Use the non linear solver to find a solution
    ceres::Solve(m_options, &m_problem, &m_summary);

    //  Print solver outputs
    std::cout << m_summary.FullReport() << "\n";
}

bool Initializer::operator()(const double* t_guess, double* t_residual) const
{
    Eigen::Matrix<double, ::CROSP::polynomial_representation::default_number_of_modes, 1> qe;
    for (unsigned int i=0;i<qe.size();i++)
        qe[i] = t_guess[i];


    Eigen::VectorXd dot_qe, ddot_qe;
    dot_qe = 0 * qe;
    ddot_qe = 0 * qe;
    m_rod_ptr->updateParameterisation(qe, dot_qe, ddot_qe);

    m_rod_ptr->forwardKinematics();
    m_rod_ptr->backwardDynamics(m_F1.block<3, 1>(0, 0),
                                m_F1.block<3, 1>(3, 0));
    const auto residual = m_rod_ptr->getStaticEquilibrium(qe);


    for (unsigned int i=0;i<residual.size();i++)
        t_residual[i] = residual[i];

    return true;
}



}   //  namespace CROSP::static_initialiser
