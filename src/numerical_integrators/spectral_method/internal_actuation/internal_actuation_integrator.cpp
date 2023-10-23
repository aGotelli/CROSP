/**
 * \file internal_actuation_integrator.cpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod internal actuator
 * \date 04-07-2023
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#include "CROSP/numerical_integrators/spectral_method/internal_actuation/internal_actuation_integrator.hpp"


namespace CROSP::numerical_integrators::spectral_method::internal_actuation_integrator {




InternalActuationIntegrator::InternalActuationIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                                         ::ATORS::tendon_driven_actuation::TendonDrivenActuation t_distributed_actuation,
                                                         const double &t_rod_length)
    : OSNI::ODEb(t_strain_parameterisation_stack->m_B_Phi_stack.cols(),
                 ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                 t_strain_parameterisation_stack->m_number_of_points,
                 t_rod_length),
      m_K_stack(t_strain_parameterisation_stack->m_K_stack),
      m_Gamma_stack(t_strain_parameterisation_stack->m_Gamma_stack),
      m_BPhi_stack(t_strain_parameterisation_stack->m_B_Phi_stack),
      m_distributed_actuation(t_distributed_actuation)
{}

InternalActuationIntegrator::InternalActuationIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack)
    : OSNI::ODEb(t_strain_parameterisation_stack->m_B_Phi_stack.cols(),
                 ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                 t_strain_parameterisation_stack->m_number_of_points)
{}


void InternalActuationIntegrator::updateActuation(const double &t_time)
{
    m_distributed_actuation.updateActuation(t_time);

    this->solveSystem();
}


Eigen::VectorXd InternalActuationIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    const Eigen::Vector3d K = m_K_stack->at(t_point);
    const Eigen::Vector3d Gamma = m_Gamma_stack->at(t_point);

    const auto cols = m_BPhi_stack.cols();
    const auto rows = 6;
    const auto start_row_index = t_point*rows;
    const auto BPhi = m_BPhi_stack.block(start_row_index, 0, rows, cols);
    const Eigen::VectorXd b_at_point =
            m_distributed_actuation.getInternalActuatedStresses(K,
                                                                Gamma,
                                                                BPhi);

    return b_at_point;

}



Eigen::VectorXd InternalActuationIntegrator::getActuation()const
{

    const Eigen::VectorXd Q_ad = this->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

    return Q_ad;
}



NullInternalActuationIntegrator::NullInternalActuationIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack)
    : InternalActuationIntegrator(t_strain_parameterisation_stack),
      m_null_actuation(Eigen::VectorXd::Zero(t_strain_parameterisation_stack->m_B_Phi_stack.cols()))
{}

}   //  namespace CROSP::numerical_integrators::spectral_method::internal_actuation_integrator

