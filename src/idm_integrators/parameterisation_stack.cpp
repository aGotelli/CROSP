#include "CROSP/idm_integrators/parameterisation_stack.hpp"

namespace CROSP {





ParameterisationStack::ParameterisationStack(const strain_parameterisation::StrainParameterisation &t_strain_parameterisation,
                                             const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                                             const unsigned int t_number_of_Chebyshev_points)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
      m_strain_parameterisation(t_strain_parameterisation),
      m_polynomial_representation(t_polynomial_representation)
{
    //  Default initialisation with nominal strain
    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension());
    updateStacks(zeros, zeros, zeros);
}




void ParameterisationStack::updateStacks(const Eigen::VectorXd &t_qe,
                                         const Eigen::VectorXd &t_dot_qe,
                                         const Eigen::VectorXd &t_ddot_qe)
{

    ::LieAlgebra::Vector6d xi;
    ::LieAlgebra::Vector6d dot_xi;
    ::LieAlgebra::Vector6d ddot_xi;


    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){

        xi      = m_map_to_strain_stack[i]*t_qe + m_strain_parameterisation.m_constant_strain;
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





}
