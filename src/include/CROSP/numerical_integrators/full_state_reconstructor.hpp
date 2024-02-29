#ifndef FULL_STATE_RECONSTRUCTOR_HPP
#define FULL_STATE_RECONSTRUCTOR_HPP

#include <iostream>
#include <Eigen/Dense>



namespace CROSP::numerical_integrators {

struct FullODEStatesObservations {

    FullODEStatesObservations()=default;

    FullODEStatesObservations(const unsigned int t_number_of_cols,
                              const unsigned int t_generalised_coordinates_dimension,
                              const std::pair<unsigned int, unsigned int> t_orientation_dimensions={4, 1});



    void writeState(const Eigen::VectorXd t_state, const unsigned int t_col);

    std::pair<unsigned int, unsigned int> m_orientation_dimensions;

    unsigned int m_generalised_coordinates_dimension;

    Eigen::MatrixXd orientation_stack;
    Eigen::MatrixXd r_stack;

    Eigen::MatrixXd Omega_stack;
    Eigen::MatrixXd V_stack;

    Eigen::MatrixXd dot_Omega_stack;
    Eigen::MatrixXd dot_V_stack;

    Eigen::MatrixXd C_stack;
    Eigen::MatrixXd N_stack;

    Eigen::MatrixXd Qa_stack;
    Eigen::MatrixXd Qad_stack;

};



}   //  namespace CROSP::numerical_integrators



#endif // FULL_STATE_RECONSTRUCTOR_HPP
