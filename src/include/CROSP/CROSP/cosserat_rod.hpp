/**
 * \file cosserat_rod.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the definition of strain based parameterised Cosserat rod
 * \date 12-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */



#ifndef COSSERAT_ROD_HPP
#define COSSERAT_ROD_HPP

#include <Eigen/Dense>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/material_properties/material_properties.hpp"

#include "CROSP/idm_integrators/idm_integrators.hpp"
#include "CROSP/tidm_integrators/tidm_integrators.hpp"

namespace CROSP {




class CosseratRod
{
public:

    CosseratRod(unsigned int t_number_of_Chebyshev_points);

    CosseratRod(std::shared_ptr<material_properties::MaterialProperties> t_material_properties,
                std::shared_ptr<strain_parameterisation::StrainParameterisation> t_strain_parameterisation,
                std::shared_ptr<strain_parameterisation::StrainParameterisationDelta> t_strain_parameterisation_perturbation);

    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                const Eigen::VectorXd &t_dot_qe,
                                const Eigen::VectorXd &t_ddot_qe);

    void forwardKinematics();

    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                           const Eigen::Vector3d &t_initial_position,
                           const Eigen::Vector3d &t_initial_angular_velocity,
                           const Eigen::Vector3d &t_initial_linear_velocity,
                           const Eigen::Vector3d &t_initial_angular_acceleration,
                           const Eigen::Vector3d &t_initial_linear_acceleration);


    void updateParameterisationPerturbation(const Eigen::VectorXd &t_Delta_qe,
                                            const Eigen::VectorXd &t_Delta_dot_qe,
                                            const Eigen::VectorXd &t_Delta_ddot_qe);

    void forwardTangentKinematics();

    void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                  const Eigen::Vector3d &t_initial_Delta_position,
                                  const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                  const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                  const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                  const Eigen::Vector3d &t_initial_Delta_linear_acceleration);

    ::LieAlgebra::Kinematics getKinematicsAtTip()const;



    void backwardDynamics(const Eigen::Vector3d &t_force_at_tip,
                          const Eigen::Vector3d &t_couple_at_tip);






    void backwardTangentDynamics(const Eigen::Vector3d &t_Delta_force_at_tip,
                                 const Eigen::Vector3d &t_Delta_couple_at_tip);




#ifndef DEVELOPER
private:
#endif

    /*!
     * \brief setForwardIntegratorsInitialConditions sets the initial conditions for the forward integrators.
     *
     * This function sets the initial conditions for the forward integrators. For these integrators, we integrate from a fixed
     * reference configuration so we can set their ivp once and for all.
     */
    void setForwardIntegratorsInitialConditions();

    std::shared_ptr<material_properties::MaterialProperties> m_material_properties { std::make_shared<material_properties::MaterialProperties>() };

    std::shared_ptr<strain_parameterisation::StrainParameterisation> m_strain_parameterisation { std::make_shared<strain_parameterisation::StrainParameterisation>() };

    std::shared_ptr<strain_parameterisation::StrainParameterisationDelta> m_strain_parameterisation_perturbation { std::make_shared<strain_parameterisation::StrainParameterisationDelta>() };

    std::shared_ptr<idm_integrators::IDMIntegrators> m_idm_integrators { std::make_shared<idm_integrators::IDMIntegrators>(m_strain_parameterisation,
                                                                                                                           m_material_properties )};

    std::shared_ptr<tidm_integrators::TIDMIntegrators> m_tidm_integrators { std::make_shared<tidm_integrators::TIDMIntegrators>(m_strain_parameterisation,
                                                                                                                                m_strain_parameterisation_perturbation,
                                                                                                                                m_idm_integrators,
                                                                                                                                m_material_properties) };
};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
