/**
 * \file tidm_integrators.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod TIDM
 * \date 28-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#ifndef TIDM_INTEGRATORS_HPP
#define TIDM_INTEGRATORS_HPP

#include <memory>
#include <vector>

#include "OSNI/include/OSNI.hpp"

#include "strain_parameterisation.hpp"
#include "idm_integrators.hpp"

namespace CROSP {



struct DeltaRotation : public OSNI::ODEAb {
    DeltaRotation(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                  std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_K_stack,
                  const unsigned int t_number_of_Chebyshev_points) :    OSNI::ODEAb(4, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                                                                                    t_number_of_Chebyshev_points),
                                                                        m_K_stack(t_K_stack),
                                                                        m_Delta_K_stack(t_Delta_K_stack)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {
        return m_Delta_K_stack->at(t_point);
    }

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
};


struct DeltaPosition : public OSNI::ODEAb {

    DeltaPosition(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                  std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                  std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_Lambda_stack,
                  std::shared_ptr<const DeltaRotation> t_delta_rotation,
                  const unsigned int t_number_of_Chebyshev_points) :    OSNI::ODEAb(4, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                                                                                    t_number_of_Chebyshev_points),
                                                                        m_K_stack(t_K_stack),
                                                                        m_Lambda_stack(t_Lambda_stack),
                                                                        m_Delta_Lambda_stack(t_Delta_Lambda_stack),
                                                                        m_delta_rotation(t_delta_rotation)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {
        return -::LieAlgebra::skew( m_Lambda_stack->at(t_point) ) * m_delta_rotation->getStateAtPoint(t_point)
                + m_Lambda_stack->at(t_point);
    }

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack;

    const std::shared_ptr<const DeltaRotation> m_delta_rotation;
};


struct TIDMIntegrators{

    TIDMIntegrators(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                    std::shared_ptr<const StrainParameterisationPerturbation> t_strain_parameterisation_perturbation,
                    std::shared_ptr<const IDMIntegrators> t_idm_integrators) :  m_strain_parameterisation(t_strain_parameterisation),
                                                                                m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation),
                                                                                m_idm_integrators(t_idm_integrators)
    {}

    TIDMIntegrators(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                    std::shared_ptr<const StrainParameterisationPerturbation> t_strain_parameterisation_perturbation,
                    std::shared_ptr<const IDMIntegrators> t_idm_integrators,
                    const unsigned int t_number_of_Chebyshev_points) :  m_strain_parameterisation(t_strain_parameterisation),
                                                                        m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation),
                                                                        m_idm_integrators(t_idm_integrators),
                                                                        m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
    {}



    const std::shared_ptr<const StrainParameterisation> m_strain_parameterisation;
    const std::shared_ptr<const StrainParameterisationPerturbation> m_strain_parameterisation_perturbation;
    const std::shared_ptr<const IDMIntegrators> m_idm_integrators;

    const unsigned int m_number_of_Chebyshev_points { 17 };


    std::shared_ptr<DeltaRotation> m_delta_rotation { std::make_shared<DeltaRotation>( m_strain_parameterisation->m_K_stack,
                                                                                       m_strain_parameterisation_perturbation->m_Delta_K_stack,
                                                                                       m_number_of_Chebyshev_points) };

    std::shared_ptr<DeltaPosition> m_delta_position { std::make_shared<DeltaPosition>( m_strain_parameterisation->m_K_stack,
                                                                                       m_strain_parameterisation->m_Lambda_stack,
                                                                                       m_strain_parameterisation_perturbation->m_Delta_Lambda_stack,
                                                                                       m_delta_rotation,
                                                                                       m_number_of_Chebyshev_points) };

};

}   //  namespace CROSP

#endif // TIDM_INTEGRATORS_HPP
