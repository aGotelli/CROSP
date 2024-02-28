/**
 * \file integration_polynomial_base.cpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod polynomial base
 * \date 04-07-2023
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#include "CROSP/numerical_integrators/spectral_method/integration_polynomial_base/integration_polynomial_base.hpp"


namespace CROSP::numerical_integrators::spectral_method::integration_polynomial_base {


unsigned int getStateSize(const Eigen::MatrixXd &t_Phi)
{
    Eigen::MatrixXd PhiTPhi = t_Phi.transpose()*t_Phi;

    return PhiTPhi.size();
}


IntegratorPolynomialBase::IntegratorPolynomialBase(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack)
    : OSNI::ODEb( getStateSize(t_strain_parameterisation_stack->m_Phi_stack[0]),
                 ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                 t_strain_parameterisation_stack->m_number_of_points),
      m_strain_parameterisation_stack(t_strain_parameterisation_stack)
{}





Eigen::VectorXd IntegratorPolynomialBase::computerParametersVectorAtPoint(const unsigned int t_point)
{

    const auto Phi = m_strain_parameterisation_stack->m_Phi_stack.at( t_point );

    Eigen::MatrixXd PhiTPhi = Phi.transpose()*Phi;

    Eigen::VectorXd b_at_point = Eigen::Map<Eigen::VectorXd>(PhiTPhi.data(), PhiTPhi.size());

    return b_at_point;

}



}   //  namespace CROSP::numerical_integrators::spectral_method::integration_polynomial_base

