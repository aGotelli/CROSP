#include "CROSP/numerical_integrators/full_state_reconstructor.hpp"



namespace CROSP::numerical_integrators {



FullODEStatesObservations::FullODEStatesObservations(const unsigned int t_number_of_cols,
                                                     const unsigned int t_generalised_coordinates_dimension,
                                                     const std::pair<unsigned int, unsigned int> t_orientation_dimensions)
    : m_orientation_dimensions(t_orientation_dimensions),
      m_generalised_coordinates_dimension(t_generalised_coordinates_dimension)
{
    orientation_stack.resize(m_orientation_dimensions.first, m_orientation_dimensions.second*t_number_of_cols);
    r_stack.resize(3, t_number_of_cols);

    Omega_stack.resize(3, t_number_of_cols);
    V_stack.resize(3, t_number_of_cols);

    dot_Omega_stack.resize(3, t_number_of_cols);
    dot_V_stack.resize(3, t_number_of_cols);

    C_stack.resize(3, t_number_of_cols);
    N_stack.resize(3, t_number_of_cols);

    Qa_stack.resize(m_generalised_coordinates_dimension, t_number_of_cols);
    Qad_stack.resize(m_generalised_coordinates_dimension, t_number_of_cols);
}



void FullODEStatesObservations::writeState(const Eigen::VectorXd t_state, const unsigned int t_col)
{

    unsigned int start_row = 0;
    const unsigned int start_col = 0;
    unsigned int block_rows = m_orientation_dimensions.first;
    const unsigned int block_cols = 1;

    orientation_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);

    //  Update row initial
    start_row += block_rows;

    block_rows = 3;

    r_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);
    //  Update row initial
    start_row += block_rows;


    Omega_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);
    //  Update row initial
    start_row += block_rows;


    V_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);
    //  Update row initial
    start_row += block_rows;


    dot_Omega_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);
    //  Update row initial
    start_row += block_rows;


    dot_V_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);
    //  Update row initial
    start_row += block_rows;


    C_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);
    //  Update row initial
    start_row += block_rows;


    N_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);
    //  Update row initial
    start_row += block_rows;


    block_rows = m_generalised_coordinates_dimension;

    Qa_stack.block(0, t_col, block_rows, 1) =
            t_state.block(start_row, start_col, block_rows, block_cols);

}



}   //  namespace CROSP::numerical_integrators



