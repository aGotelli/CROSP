#ifndef INTERNALLYACTUATEDCOSSERATROD_HPP
#define INTERNALLYACTUATEDCOSSERATROD_HPP

#include "cosserat_rod.hpp"

namespace CROSP {


struct Tendon{

    Tendon(const Eigen::Vector3d &t_D,
           const Eigen::Vector3d &t_D_prime=Eigen::Vector3d::Zero())
        : m_D(t_D),
          m_D_prime(t_D_prime)
    {}


    Eigen::Vector3d m_D;
    Eigen::Vector3d m_D_prime;

    ::LieAlgebra::Vector6d getStrain(const double &t_tau,
                                     const Eigen::Vector3d &t_K,
                                     const Eigen::Vector3d &t_Gamma)
    {
        const auto K_skew = ::LieAlgebra::skew(t_K);
        const auto D_skew = ::LieAlgebra::skew(m_D);

        const Eigen::Vector3d T = t_Gamma + K_skew*m_D + m_D_prime;

        ::LieAlgebra::Vector6d strain_map;
        strain_map <<  D_skew*T,
                           T   ;


        const double norm = T.norm();

        const ::LieAlgebra::Vector6d strain_from_cable = strain_map*t_tau / norm;

        return strain_from_cable;
    }


};

struct Actuator {
    std::function<void(const double&)> updateTorque;

    double m_tau;
};

typedef std::pair<Actuator, Tendon> ActuatedTendon;


struct InternalActuation{


    void addTendon(const Tendon &t_tendon,
                   const Actuator &t_actuator);


    std::vector<ActuatedTendon> m_actuated_tendons;

};

struct GeneralisedInternalActuationIntegrator : public OSNI::ODEb {

    GeneralisedInternalActuationIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                           InternalActuation t_internal_actuation);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    std::shared_ptr<const strain_parameterisation::StrainParameterisation> m_strain_parameterisation;

    InternalActuation m_internal_actuation;

    Eigen::VectorXd m_tau;

};



class InternallyActuatedCosseratRod : public CosseratRod
{
public:
    InternallyActuatedCosseratRod();


    /*!
     * \brief InternallyActuatedCosseratRod initialised the rod by giving its properties anly
     * \param t_rod_properties is the sef of desiderd properties as a ::CROSP::rod_properties::RodProperties object
     *
     * In this constructo the strain parameterisation is initialised by its defaul value
     */
    InternallyActuatedCosseratRod(const std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                  const InternalActuation t_actuation);


    /*!
     * \brief InternallyActuatedCosseratRod initialised the rod by giving its strain parameterisation
     * \param t_strain_parameterisation is the sef of desiderd parameterisation of the rod strain field as a ::CROSP::strain_parameterisation::StrainParameterisation object
     *
     * In this constructo the rod properties is initialised by its defaul value
     */
    InternallyActuatedCosseratRod(const std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  const InternalActuation t_actuation);


    /*!
     * \brief InternallyActuatedCosseratRod this constructor requires both argument to initialise a custom Cosserat rod
     * \param t_strain_parameterisation is the sef of desiderd parameterisation of the rod strain field as a ::CROSP::strain_parameterisation::StrainParameterisation object
     * \param t_rod_properties is the sef of desiderd properties as a ::CROSP::rod_properties::RodProperties object
     */
    InternallyActuatedCosseratRod(const std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                                  const std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                                  const InternalActuation t_actuation);


    virtual void updateInternalActuation(const double &t_current_time) override;

    virtual Eigen::VectorXd getQad() override;


    GeneralisedInternalActuationIntegrator m_internal_actuation_integrator;


};

}   //  namespace CROSP

#endif // INTERNALLYACTUATEDCOSSERATROD_HPP
