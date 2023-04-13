#ifndef INTERNALLYACTUATEDCOSSERATROD_HPP
#define INTERNALLYACTUATEDCOSSERATROD_HPP

#include "cosserat_rod.hpp"

#include "ATORS/tendon_driven_actuation/tendon_driven_actuation.hpp"

namespace CROSP {



class InternallyActuatedCosseratRod : public CosseratRod
{
public:
    InternallyActuatedCosseratRod();


    /*!
     * \brief InternallyActuatedCosseratRod initialised the rod by giving its properties anly
     * \param t_actuated_tendons is the set of actuated tendons acting on the rod
     *
     * In this constructo the strain parameterisation is initialised by its defaul value
     */
    InternallyActuatedCosseratRod(const ::ATORS::tendon_driven_actuation::ActuatedTendons t_actuated_tendons);



    /*!
     * \brief InternallyActuatedCosseratRod initialised the rod by giving its properties anly
     * \param t_rod_properties is the sef of desiderd properties as a ::CROSP::rod_properties::RodProperties object
     * \param t_actuated_tendons is the set of actuated tendons acting on the rod
     *
     * In this constructo the strain parameterisation is initialised by its defaul value
     */
    InternallyActuatedCosseratRod(const std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                  const ::ATORS::tendon_driven_actuation::ActuatedTendons t_actuated_tendons);


    /*!
     * \brief InternallyActuatedCosseratRod initialised the rod by giving its strain parameterisation
     * \param t_strain_parameterisation is the sef of desiderd parameterisation of the rod strain field as a ::CROSP::strain_parameterisation::StrainParameterisation object
     * \param t_actuated_tendons is the set of actuated tendons acting on the rod
     *
     * In this constructo the rod properties is initialised by its defaul value
     */
    InternallyActuatedCosseratRod(const std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  const ::ATORS::tendon_driven_actuation::ActuatedTendons t_actuated_tendons);


    /*!
     * \brief InternallyActuatedCosseratRod this constructor requires both argument to initialise a custom Cosserat rod
     * \param t_strain_parameterisation is the sef of desiderd parameterisation of the rod strain field as a ::CROSP::strain_parameterisation::StrainParameterisation object
     * \param t_rod_properties is the sef of desiderd properties as a ::CROSP::rod_properties::RodProperties object
     * \param t_actuated_tendons is the set of actuated tendons acting on the rod
     */
    InternallyActuatedCosseratRod(const std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  const std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                  const ::ATORS::tendon_driven_actuation::ActuatedTendons t_actuated_tendons);


    virtual void updateInternalActuation(const double &t_current_time) override;

    virtual Eigen::VectorXd getQad() const override;


    ::ATORS::tendon_driven_actuation::TendonDrivenActuation m_tendon_driven_actuator;


};

}   //  namespace CROSP

#endif // INTERNALLYACTUATEDCOSSERATROD_HPP
