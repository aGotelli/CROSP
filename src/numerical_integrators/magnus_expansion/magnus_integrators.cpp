#include "CROSP/numerical_integrators/magnus_expansion/magnus_integrators.hpp"


namespace CROSP::numerical_integrators::magnus_expansion {


MagnusIntegrators::MagnusIntegrators(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                  std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_Delta_strain_parameterisation_stack,
                  std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
    : m_strain_parameterisation_stack(t_strain_parameterisation_stack),
      m_Delta_strain_parameterisation_stack(t_Delta_strain_parameterisation_stack),
      m_rod_properties(t_rod_properties)
{}

}   //  namespace CROSP::numerical_integrators::magnus_expansion
