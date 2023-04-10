#include "CROSP/CROSP/internally_actuated_cosserat_rod.hpp"



namespace CROSP {



void InternalActuation::addTendon(const Tendon &t_tendon,
                                  const Actuator &t_actuator)
{
    m_actuated_tendons.push_back( {t_actuator, t_tendon} );
}



GeneralisedInternalActuationIntegrator::GeneralisedInternalActuationIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                                               InternalActuation t_internal_actuation)
    : OSNI::ODEb(t_strain_parameterisation->m_polynomial_representation->getCoordinatesDimension(),
                 ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                 t_strain_parameterisation->m_number_of_Chebyshev_points),
      m_strain_parameterisation(t_strain_parameterisation),
      m_internal_actuation(t_internal_actuation)
{}

Eigen::VectorXd GeneralisedInternalActuationIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    const Eigen::Vector3d K = m_strain_parameterisation->m_K_stack->at(t_point);
    const Eigen::Vector3d Gamma = m_strain_parameterisation->m_Gamma_stack->at(t_point);


    ::LieAlgebra::Vector6d strain_from_cables = ::LieAlgebra::Vector6d::Zero();

    for(auto& [actuator, tendon] : m_internal_actuation.m_actuated_tendons){
        const auto tau = actuator.m_tau;

        strain_from_cables += tendon.getStrain(tau, K, Gamma);
    }


    const auto BPhi_T = m_strain_parameterisation->m_map_to_strain_stack[t_point].transpose();
    const Eigen::VectorXd b_at_point = BPhi_T * strain_from_cables;

    return b_at_point;

}




InternallyActuatedCosseratRod::InternallyActuatedCosseratRod(const std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                             const InternalActuation t_actuation)
    : CosseratRod(t_rod_properties),
      m_internal_actuation_integrator(this->m_strain_parameterisation,
                                      t_actuation)
{}



InternallyActuatedCosseratRod::InternallyActuatedCosseratRod(const std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                             const InternalActuation t_actuation)
    : CosseratRod(t_strain_parameterisation),
      m_internal_actuation_integrator(this->m_strain_parameterisation,
                                      t_actuation)
{}



InternallyActuatedCosseratRod::InternallyActuatedCosseratRod(const std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                                             const std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                                             const InternalActuation t_actuation)
    : CosseratRod(t_strain_parameterisation, t_rod_properties),
      m_internal_actuation_integrator(this->m_strain_parameterisation,
                                      t_actuation)
{}





void InternallyActuatedCosseratRod::updateInternalActuation(const double &t_current_time)
{
    std::for_each(m_internal_actuation_integrator.m_internal_actuation.m_actuated_tendons.begin(),
                  m_internal_actuation_integrator.m_internal_actuation.m_actuated_tendons.end(),
                  [&](ActuatedTendon &t_actuated_tendon){
        t_actuated_tendon.first.updateTorque(t_current_time);
    });
}

Eigen::VectorXd InternallyActuatedCosseratRod::getQad()
{
    m_internal_actuation_integrator.solveSystem();

    const auto Q_ad = m_internal_actuation_integrator.getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

    return Q_ad;

}










}   //  namespace CROSP
