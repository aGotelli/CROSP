#include "CROSP/CROSP/internally_actuated_cosserat_rod.hpp"



namespace CROSP {


InternallyActuatedCosseratRod::InternallyActuatedCosseratRod(const ::ATORS::tendon_driven_actuation::ActuatedTendons t_actuated_tendons)
    : CosseratRod(),
      m_tendon_driven_actuator(this->m_strain_parameterisation->m_K_stack,
                               this->m_strain_parameterisation->m_Gamma_stack,
                               this->m_strain_parameterisation->m_map_to_strain_stack,
                               t_actuated_tendons,
                               m_rod_properties->m_rod_dimensions.m_L)
{}


InternallyActuatedCosseratRod::InternallyActuatedCosseratRod(const std::shared_ptr<rod_properties::RodProperties> t_rod_properties,
                                                             const ::ATORS::tendon_driven_actuation::ActuatedTendons t_actuated_tendons)
    : CosseratRod(t_rod_properties),
      m_tendon_driven_actuator(this->m_strain_parameterisation->m_K_stack,
                               this->m_strain_parameterisation->m_Gamma_stack,
                               this->m_strain_parameterisation->m_map_to_strain_stack,
                               t_actuated_tendons,
                               m_rod_properties->m_rod_dimensions.m_L)
{}



InternallyActuatedCosseratRod::InternallyActuatedCosseratRod(const std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                             const ::ATORS::tendon_driven_actuation::ActuatedTendons t_actuated_tendons)
    : CosseratRod(t_strain_parameterisation),
      m_tendon_driven_actuator(this->m_strain_parameterisation->m_K_stack,
                               this->m_strain_parameterisation->m_Gamma_stack,
                               this->m_strain_parameterisation->m_map_to_strain_stack,
                               t_actuated_tendons,
                               m_rod_properties->m_rod_dimensions.m_L)
{}



InternallyActuatedCosseratRod::InternallyActuatedCosseratRod(const std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                             const std::shared_ptr<rod_properties::RodProperties> t_rod_properties,
                                                             const ::ATORS::tendon_driven_actuation::ActuatedTendons t_actuated_tendons)
    : CosseratRod(t_strain_parameterisation, t_rod_properties),
      m_tendon_driven_actuator(this->m_strain_parameterisation->m_K_stack,
                               this->m_strain_parameterisation->m_Gamma_stack,
                               this->m_strain_parameterisation->m_map_to_strain_stack,
                               t_actuated_tendons,
                               m_rod_properties->m_rod_dimensions.m_L)
{}


void InternallyActuatedCosseratRod::updateInternalActuation(const double &t_current_time)
{
    m_tendon_driven_actuator.updateActuation(t_current_time);
}

Eigen::VectorXd InternallyActuatedCosseratRod::getQad()const
{
    return m_tendon_driven_actuator.getActuation();
}

void InternallyActuatedCosseratRod::updateRodLength(const double &t_rod_lenght)
{
    m_rod_properties->updateRodLength( t_rod_lenght );

    m_idm_integrators->updateIntegrationDomain( t_rod_lenght );

    m_tidm_integrators->updateIntegrationDomain( t_rod_lenght );

    m_tendon_driven_actuator.setUpperIntegrationDomain( t_rod_lenght );
}


}   //  namespace CROSP
