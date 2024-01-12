#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"
#include <math_tools/Chebyshev/chebyshev_differentiation.hpp>

namespace CROSP::strain_parameterisation_stack {





StrainParameterisationStack::StrainParameterisationStack(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                                         const unsigned int t_number_of_Chebyshev_point,
                                                         const StrainFunction &t_Xi_c)
    : StrainParameterisationStack(t_polynomial_representation,
                                  ::Chebyshev::ComputeChebyshevPoints(t_number_of_Chebyshev_point),
                                  t_Xi_c)
{}

StrainParameterisationStack::StrainParameterisationStack(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                                         const std::vector<double> &t_observation_points,
                                                         const StrainFunction &t_Xi_c)
    : m_number_of_points(t_observation_points.size()),
      m_Xi_c(t_Xi_c)
{

    const Eigen::MatrixXd B = t_polynomial_representation.m_B;
    Eigen::MatrixXd Phi = t_polynomial_representation.getPhi( 0 );
    Eigen::MatrixXd BPhi = B*Phi;


    m_B_Phi_stack =
        Eigen::MatrixXd::Zero(m_number_of_points*BPhi.rows(), BPhi.cols());

    m_Xi_stack =
            Eigen::MatrixXd::Zero(m_number_of_points*6, 1);

    m_dot_Xi_stack =
            Eigen::MatrixXd::Zero(m_number_of_points*6, 1);

    m_ddot_Xi_stack =
            Eigen::MatrixXd::Zero(m_number_of_points*6, 1);

    m_Xi_c_stack =
            Eigen::MatrixXd::Zero(m_number_of_points*6, 1);


    ::LieAlgebra::Vector6d Xi_c = ::LieAlgebra::Vector6d::Zero();
    Xi_c(3) = 1;

    unsigned int start_row;
    constexpr unsigned int block_row_se3 = 6;
    constexpr unsigned int block_row_vector = 3;
    for(unsigned int step=0; const auto point : t_observation_points){

            //  Compute the base
            Phi = t_polynomial_representation.getPhi( point );

            //  Add base to stack
            m_Phi_stack.push_back(Phi);

            //  Compute corresponding row
            start_row = step * block_row_se3;

            //  Add BPhi to matrix
            m_B_Phi_stack.block(start_row, 0, block_row_se3, BPhi.cols()) = B*Phi;

            //  Add all the constrained strains in matrix form
            m_Xi_c_stack.block(start_row, 0, block_row_se3, 1) = m_Xi_c(point);


            start_row = block_row_se3*step;
            m_K_stack->push_back(           m_Xi_stack.block<block_row_vector, 1>(start_row, 0) );
            m_dot_K_stack->push_back(   m_dot_Xi_stack.block<block_row_vector, 1>(start_row, 0) );
            m_ddot_K_stack->push_back( m_ddot_Xi_stack.block<block_row_vector, 1>(start_row, 0) );

            start_row += block_row_vector;
            m_Gamma_stack->push_back(           m_Xi_stack.block<block_row_vector, 1>(start_row, 0) );
            m_dot_Gamma_stack->push_back(   m_dot_Xi_stack.block<block_row_vector, 1>(start_row, 0) );
            m_ddot_Gamma_stack->push_back( m_ddot_Xi_stack.block<block_row_vector, 1>(start_row, 0) );


            step++;
    }


}


void StrainParameterisationStack::updateStrainParameterisation(const Eigen::VectorXd &t_q,
                                                               const Eigen::VectorXd &t_dot_q,
                                                               const Eigen::VectorXd &t_ddot_q)
{
    for(unsigned int step=0; step<m_number_of_points; step++){

        std::cout << m_Gamma_stack->at(step) << "\n\n";
    }


    m_Xi_stack = m_B_Phi_stack * t_q + m_Xi_c_stack;
    m_dot_Xi_stack = m_B_Phi_stack * t_dot_q;
    m_ddot_Xi_stack = m_B_Phi_stack * t_ddot_q;

    std::cout << "Update in Xi : \n" << m_Xi_stack << "\n\n\n";


    std::cout << "\n now K stack : \n";
    for(unsigned int step=0; step<m_number_of_points; step++){

        std::cout << m_Gamma_stack->at(step) << "\n\n";
    }
}





}   //  namespace CROSP::strain_parameterisation_stack
