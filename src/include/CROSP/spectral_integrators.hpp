#ifndef SPECTRAL_INTEGRATORS_HPP
#define SPECTRAL_INTEGRATORS_HPP

#include "idm_integrators/idm_integrators.hpp"
#include "tidm_integrators/tidm_integrators.hpp"

namespace CROSP::integrators::spectral_integrators {


struct SpectralIntegrators{


    SpectralIntegrators(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                        std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
                        std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
        :   m_idm_integrators( std::make_shared<idm_integrators::IDMIntegrators>(t_strain_parameterisation,
                                                                                 t_rod_properties) ),
            m_tidm_integrators( std::make_shared<tidm_integrators::TIDMIntegrators>(t_strain_parameterisation,
                                                                                    t_strain_parameterisation_Delta,
                                                                                    m_idm_integrators,
                                                                                    t_rod_properties) )
    {}




    //  The set of integrators needed for the IDM
    std::shared_ptr<idm_integrators::IDMIntegrators> m_idm_integrators;

    //  The set of integrators needed for the TIDM
    std::shared_ptr<tidm_integrators::TIDMIntegrators> m_tidm_integrators;



};


//SpectralIntegrators::SpectralIntegrators(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
//                                         std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation_Delta,
//                                         std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
//    :   m_idm_integrators( std::make_shared<idm_integrators::IDMIntegrators>(t_strain_parameterisation,
//                                                                             t_rod_properties) ),
//        m_tidm_integrators( std::make_shared<tidm_integrators::TIDMIntegrators>(t_strain_parameterisation,
//                                                                                t_strain_parameterisation_Delta,
//                                                                                m_idm_integrators,
//                                                                                t_rod_properties) )
//{}



typedef std::shared_ptr<SpectralIntegrators> IntegratorsSPtr;




}   //  namespace CROSP::integrators::spectral_integrators

#endif // SPECTRAL_INTEGRATORS_HPP
