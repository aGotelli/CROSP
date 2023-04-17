#ifndef INTERNALACTUATION_HPP
#define INTERNALACTUATION_HPP

#include "OSNI/ODE/ode_b.hpp"
#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include <memory>

namespace CROSP::internal_actuation {

struct GeneralisedInternalActuationIntegrator : public OSNI::ODEb {

    GeneralisedInternalActuationIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;

    Eigen::VectorXd m_tau;

};




class InternalActuation
{
public:
    InternalActuation(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation);

    virtual void setActuation(const double &t_current_time);

    Eigen::VectorXd computedGeneraliseInternalActuation();

    GeneralisedInternalActuationIntegrator m_integrator;
};


}   //  namespace CROSP::internal_actuation

#endif // INTERNALACTUATION_HPP
