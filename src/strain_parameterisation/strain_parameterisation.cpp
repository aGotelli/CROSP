#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace CROSP::strain_parameterisation {




StrainParameterisation::StrainParameterisation(const unsigned int t_number_of_Chebyshev_points)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{
    //  Default initialisation with nominal strain
    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(m_polynomial_representation->getCoordinatesDimension());
    updateStacks(zeros, zeros, zeros);
}


StrainParameterisation::StrainParameterisation(const std::shared_ptr<const polynomial_representation::PolynomialRepresentation> t_polynomial_representation,
                                               const unsigned int t_number_of_Chebyshev_points)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_polynomial_representation(t_polynomial_representation)
{
    //  Default initialisation with nominal strain
    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(m_polynomial_representation->getCoordinatesDimension());
    updateStacks(zeros, zeros, zeros);
}



StrainParameterisation::StrainParameterisation(const std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                                               const Eigen::VectorXd t_constant_strain)
    : m_number_of_Chebyshev_points(t_strain_parameterisation->m_number_of_Chebyshev_points),
      m_constant_strain(t_constant_strain),
      m_polynomial_representation(t_strain_parameterisation->m_polynomial_representation)
{
    //  Default initialisation with nominal strain
    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(m_polynomial_representation->getCoordinatesDimension());
    updateStacks(zeros, zeros, zeros);
}




void StrainParameterisation::updateStacks(const Eigen::VectorXd &t_qe,
                                          const Eigen::VectorXd &t_dot_qe,
                                          const Eigen::VectorXd &t_ddot_qe)
{

    ::LieAlgebra::Vector6d xi;
    ::LieAlgebra::Vector6d dot_xi;
    ::LieAlgebra::Vector6d ddot_xi;


    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){

        xi      = m_map_to_strain_stack[i]*t_qe + m_constant_strain;
        dot_xi  = m_map_to_strain_stack[i]*t_dot_qe;
        ddot_xi = m_map_to_strain_stack[i]*t_ddot_qe;

        m_K_stack->at(i)      = xi.block<3,1>(0,0);
        m_dot_K_stack->at(i)  = dot_xi.block<3,1>(0,0);
        m_ddot_K_stack->at(i) = ddot_xi.block<3,1>(0,0);

        m_Gamma_stack->at(i)      = xi.block<3,1>(3,0);
        m_dot_Gamma_stack->at(i)  = dot_xi.block<3,1>(3,0);
        m_ddot_Gamma_stack->at(i) = ddot_xi.block<3,1>(3,0);

    }

}



}   //  namespace CROSP::StrainParameterisation
