#ifndef INTERNALACTUATIONINTEGRATOR_HPP
#define INTERNALACTUATIONINTEGRATOR_HPP


#include <memory>

#include "CROSP/rod_properties/rod_properties.hpp"

#include "OSNI/OSNI.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"

#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"
#include "ATORS/ATORS.hpp"


/// \brief CROSP::idm_integrators is the namespace containing the integrators for the internal actuation
namespace CROSP::numerical_integrators::spectral_method::internal_actuation_integrator {


struct InternalActuationIntegrator : public ::ATORS::ActuatorBase, ::OSNI::ODEb{




    InternalActuationIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                                ::ATORS::distributed_actuation::DistributedActuationUptr t_distributed_actuation,
                                const double &t_rod_length);

protected:
    InternalActuationIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack);
public:

    virtual void updateActuation(const double &t_time)override;

    virtual void setActuation(const Eigen::VectorXd &)override{}

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)override;


    virtual Eigen::VectorXd getActuation()const override;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack;

    Eigen::MatrixXd m_BPhi_stack;



    ::ATORS::distributed_actuation::DistributedActuationUptr m_distributed_actuation;


};

struct NullInternalActuationIntegrator : public InternalActuationIntegrator{




    NullInternalActuationIntegrator(std::shared_ptr<const strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack);

    virtual void updateActuation(const double &)final{};

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int )final{ return m_null_actuation;};

    virtual void setActuation(const Eigen::VectorXd &)final{}

    virtual Eigen::VectorXd getActuation()const final{ return m_null_actuation; };


    Eigen::VectorXd m_null_actuation;

};


typedef std::unique_ptr<InternalActuationIntegrator> InternalActuationIntegratorUptr;

typedef std::unique_ptr<NullInternalActuationIntegrator> NullInternalActuationIntegratorUptr;



}   //  namespace CROSP::numerical_integrators::spectral_method::internal_actuation_integrator

#endif // INTERNALACTUATIONINTEGRATOR_HPP
