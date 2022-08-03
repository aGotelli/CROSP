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

#include "strain_parameterisation.hpp"
#include "material_properties.hpp"

#include "idm_integrators.hpp"
#include "tidm_integrators.hpp"

namespace CROSP {





class CosseratRod
{
public:
    CosseratRod()=default;

    CosseratRod(unsigned int t_number_of_Chebyshev_points);

    CosseratRod(std::shared_ptr<MaterialProperties> t_material_properties,
                std::shared_ptr<StrainParameterisation> t_strain_parameterisation,
                std::shared_ptr<StrainParameterisationPerturbation> t_strain_parameterisation_perturbation) :   m_material_properties(t_material_properties),
                                                                                                                m_strain_parameterisation(t_strain_parameterisation),
                                                                                                                m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation)
    {}

    void updateParameterisation(const Eigen::VectorXd &t_qe,
                                const Eigen::VectorXd &t_dot_qe,
                                const Eigen::VectorXd &t_ddot_qe);

    void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion=Eigen::Vector4d(1, 0, 0, 0),
                           const Eigen::Vector3d &t_initial_position=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_angular_velocity=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_linear_velocity=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_angular_acceleration=Eigen::Vector3d::Zero(),
                           const Eigen::Vector3d &t_initial_linear_acceleration=Eigen::Vector3d::Zero());


    void updateParameterisationPerturbation(const Eigen::VectorXd &t_Delta_qe,
                                            const Eigen::VectorXd &t_Delta_dot_qe,
                                            const Eigen::VectorXd &t_Delta_ddot_qe);

    void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation=Eigen::Vector3d::Zero(),
                                  const Eigen::Vector3d &t_initial_Delta_position=Eigen::Vector3d::Zero(),
                                  const Eigen::Vector3d &t_initial_Delta_angular_velocity=Eigen::Vector3d::Zero(),
                                  const Eigen::Vector3d &t_initial_Delta_linear_velocity=Eigen::Vector3d::Zero(),
                                  const Eigen::Vector3d &t_initial_Delta_angular_acceleration=Eigen::Vector3d::Zero(),
                                  const Eigen::Vector3d &t_initial_Delta_linear_acceleration=Eigen::Vector3d::Zero());

    ::LieAlgebra::Kinematics getKinematicsAtTip()const;



    void backwardDynamics(const Eigen::Vector3d &t_force_at_tip,
                          const Eigen::Vector3d &t_couple_at_tip);



#ifndef DEVELOPER
private:
#endif

    std::shared_ptr<MaterialProperties> m_material_properties { std::make_shared<MaterialProperties>() };

    std::shared_ptr<StrainParameterisation> m_strain_parameterisation { std::make_shared<StrainParameterisation>() };

    std::shared_ptr<StrainParameterisationPerturbation> m_strain_parameterisation_perturbation { std::make_shared<StrainParameterisationPerturbation>() };

    std::shared_ptr<IDMIntegrators> m_idm_integrators { std::make_shared<IDMIntegrators>(m_strain_parameterisation,
                                                                                         m_material_properties )};

    std::shared_ptr<TIDMIntegrators> m_tidm_integrators { std::make_shared<TIDMIntegrators>(m_strain_parameterisation,
                                                                                            m_strain_parameterisation_perturbation,
                                                                                            m_idm_integrators,
                                                                                            m_material_properties) };
};



}   //  namespace CROSP

#endif // COSSERAT_ROD_HPP
