/**
 * \file polynomial_representation.cpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the definition of class and function used for the parameterisation of the strain field with a polynomial basis
 * \date 12-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */

#include "CROSP/polynomial_representation/polynomial_representation.hpp"


namespace CROSP::polynomial_representation {




PolynomialRepresentation::PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations,
                                                   const unsigned int t_number_of_modes,
                                                   const PolynomialBase t_polynomial_base)
    : m_admitted_deformations(t_admitted_deformations),
      m_number_of_modes_stack( [&](){
                    std::vector<unsigned int> number_of_modes_stack;

                    //  Populate the vector for every true entry of the admitted deformations
                    std::for_each(t_admitted_deformations.begin(),
                                  t_admitted_deformations.end(),
                                  [&](const bool dof){
                        if(dof == true)
                            number_of_modes_stack.push_back( t_number_of_modes );
                    });
                    return number_of_modes_stack;}() ),
      m_polynomial_base(t_polynomial_base)
{}




PolynomialRepresentation::PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations,
                                                   const std::vector<unsigned int> t_number_of_modes_stack,
                                                   const PolynomialBase t_polynomial_base)
    : m_admitted_deformations(t_admitted_deformations),
      m_number_of_modes_stack(t_number_of_modes_stack),
      m_polynomial_base(t_polynomial_base)
{
    if(m_number_of_modes_stack.size() != m_na){
        std::string help_message = "The array containing the number of modes must have as many"
                                   " elements as the allowed deformations. Here you passed : " + std::to_string( m_number_of_modes_stack.size() ) +
                                   " entries for the number of modes stack while you activated : " + std::to_string( m_na ) +
                                   " deformations";


            throw std::runtime_error{help_message};
    }
}


PolynomialRepresentation::PolynomialRepresentation(const PolynomialRepresentation &t_other)
    : m_admitted_deformations(t_other.m_admitted_deformations),
      m_number_of_modes_stack(t_other.m_number_of_modes_stack),
      m_polynomial_base(t_other.m_polynomial_base)
{}




Eigen::MatrixXd PolynomialRepresentation::getPhi(const double& t_X,
                                                 const double& t_begin,
                                                 const double& t_end) const
{
    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    double x = ( 2 * t_X - ( t_end + t_begin) ) / ( t_end - t_begin );


    //  Define the Phi matrix beforehand
    Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(m_na, m_total_number_of_modes);

    for(unsigned int col_index=0, row_index = 0; const auto &number_of_modes : m_number_of_modes_stack){

        //  Create the vector
        Eigen::VectorXd Phi_i(number_of_modes);
        for(unsigned int i=0; i<number_of_modes; i++)
            Phi_i[i] = m_polynomial_base(i, x);

        //  Put it in the matrix
        Phi.block(row_index, col_index, 1, number_of_modes) = Phi_i.transpose();

        row_index ++;
        col_index += number_of_modes;
    }


    return Phi;
}








}   //  namespace CROSP::polynomial_representation
