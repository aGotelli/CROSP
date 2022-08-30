#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace CROSP::strain_parameterisation {




StrainParameterisation::StrainParameterisation(const unsigned int t_number_of_Chebyshev_points)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


StrainParameterisation::StrainParameterisation(const Eigen::VectorXd &t_constant_strain)
    : m_constrained_strain( defineConstrainedStrain(t_constant_strain) )
{}


StrainParameterisation::StrainParameterisation(const polynomial_representation::PolynomialRepresentation t_polynomial_representation)
    : m_polynomial_representation(t_polynomial_representation)
{}


StrainParameterisation::StrainParameterisation(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                               const unsigned int t_number_of_Chebyshev_points)
    : m_polynomial_representation(t_polynomial_representation),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


StrainParameterisation::StrainParameterisation(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                               const Eigen::VectorXd &t_constant_strain)
    : m_polynomial_representation(t_polynomial_representation),
      m_constrained_strain( defineConstrainedStrain(t_constant_strain) )
{}


StrainParameterisation::StrainParameterisation(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                               const Eigen::VectorXd &t_constant_strain,
                                               const unsigned int t_number_of_Chebyshev_points)
    : m_polynomial_representation(t_polynomial_representation),
      m_constrained_strain( defineConstrainedStrain(t_constant_strain) ),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}



void StrainParameterisation::updateStacks(const Eigen::VectorXd &t_qe,
                                          const Eigen::VectorXd &t_dot_qe,
                                          const Eigen::VectorXd &t_ddot_qe)
{

    Eigen::VectorXd xi;
    Eigen::VectorXd dot_xi;
    Eigen::VectorXd ddot_xi;

    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){

        xi = m_map_to_strain_stack[i]*t_qe + m_polynomial_representation.m_Bbar*m_constrained_strain;
        dot_xi = m_map_to_strain_stack[i]*t_dot_qe;
        ddot_xi = m_map_to_strain_stack[i]*t_ddot_qe;

        m_K_stack->at(i) = xi.block<3,1>(0,0);
        m_dot_K_stack->at(i) = dot_xi.block<3,1>(0,0);
        m_ddot_K_stack->at(i) = ddot_xi.block<3,1>(0,0);

        m_Lambda_stack->at(i) = xi.block<3,1>(3,0);
        m_dot_Lambda_stack->at(i) = dot_xi.block<3,1>(3,0);
        m_ddot_Lambda_stack->at(i) = ddot_xi.block<3,1>(3,0);

    }

}



Eigen::VectorXd StrainParameterisation::defineConstrainedStrain(const Eigen::VectorXd &t_constant_strain)const
{
    std::vector<int> indexes;
    std::for_each(m_polynomial_representation.m_admitted_deformations.begin(),
                  m_polynomial_representation.m_admitted_deformations.end(),
                  [&indexes, index=0](const bool dof)mutable{   if(dof == false)
                                                                    indexes.push_back(index);
                                                                index++;});
    Eigen::VectorXd xi_c = t_constant_strain(indexes);

    return xi_c;
}








}   //  namespace CROSP::StrainParameterisation
