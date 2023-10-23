#ifndef INTEGRATION_POLYNOMIAL_BASE_HPP
#define INTEGRATION_POLYNOMIAL_BASE_HPP


#include <memory>

#include "CROSP/rod_properties/rod_properties.hpp"
#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include "OSNI/OSNI.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"

#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"



/// \brief CROSP::idm_integrators is the namespace containing the integrators for the internal actuation
namespace CROSP::numerical_integrators::spectral_method::integration_polynomial_base {


struct IntegratorPolynomialBase : public ::OSNI::ODEb{


    IntegratorPolynomialBase(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack);


    std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)override;


};






}   //  namespace CROSP::numerical_integrators::spectral_method::integration_polynomial_base

#endif // INTEGRATION_POLYNOMIAL_BASE_HPP
