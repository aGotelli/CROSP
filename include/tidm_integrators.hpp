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

#include "OSNI/OSNI.hpp"

#include "material_properties.hpp"

#include "strain_parameterisation.hpp"
#include "idm_integrators.hpp"

namespace CROSP {



struct DeltaRotation : public OSNI::ODEAb {
    DeltaRotation(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                  std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_K_stack,
                  const unsigned int t_number_of_Chebyshev_points) :    OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
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
                  std::shared_ptr<const OSNI::ODESolverInterface> t_delta_rotation,
                  const unsigned int t_number_of_Chebyshev_points) :    OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
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
                + m_Delta_Lambda_stack->at(t_point);
    }

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_delta_rotation;
};




struct DeltaAngularVelocity : public OSNI::ODEAb {
    DeltaAngularVelocity(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                         std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_K_stack,
                         std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_dot_K_stack,
                         std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                         const unsigned int t_number_of_Chebyshev_points) :     OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                                                                                            t_number_of_Chebyshev_points),
                                                                                m_K_stack(t_K_stack),
                                                                                m_Delta_K_stack(t_Delta_K_stack),
                                                                                m_Delta_dot_K_stack(t_Delta_dot_K_stack),
                                                                                m_angular_velocity_integrator(t_angular_velocity_integrator)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {
        return -::LieAlgebra::skew( m_Delta_K_stack->at(t_point) )*m_angular_velocity_integrator->getStateAtPoint(t_point)
         +m_Delta_dot_K_stack->at(t_point);
    }

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;

};


struct DeltaLinearVelocity : public OSNI::ODEAb {
    DeltaLinearVelocity(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                        std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_K_stack,
                        std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Lambda_stack,
                        std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_Lambda_stack,
                        std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_dot_Lambda_stack,
                        std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                        std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                        std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                        const unsigned int t_number_of_Chebyshev_points) :  OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                                                                                        t_number_of_Chebyshev_points),
                                                                            m_K_stack(t_K_stack),
                                                                            m_Delta_K_stack(t_Delta_K_stack),
                                                                            m_Lambda_stack(t_Lambda_stack),
                                                                            m_Delta_Lambda_stack(t_Delta_Lambda_stack),
                                                                            m_Delta_dot_Lambda_stack(t_Delta_dot_Lambda_stack),
                                                                            m_angular_velocity_integrator(t_angular_velocity_integrator),
                                                                            m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
                                                                            m_linear_velocity_integrator(t_linear_velocity_integrator)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {
        return m_Delta_dot_Lambda_stack->at(t_point)
                - ::LieAlgebra::skew( m_Lambda_stack->at(t_point) ) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ) * m_linear_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_Lambda_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point);
    }



    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_Lambda_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;

};


struct DeltaAngularAccelerations : public OSNI::ODEAb {

    DeltaAngularAccelerations(std::shared_ptr<const std::vector<Eigen::Vector3d>> t_K_stack,
                              std::shared_ptr<const std::vector<Eigen::Vector3d>> t_dot_K_stack,
                              std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_K_stack,
                              std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_dot_K_stack,
                              std::shared_ptr<const std::vector<Eigen::Vector3d>> t_Delta_ddot_K_stack,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                              const unsigned int t_number_of_Chebyshev_points) :  OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                                                                                        t_number_of_Chebyshev_points),
                                                                                  m_K_stack(t_K_stack),
                                                                                  m_dot_K_stack(t_dot_K_stack),
                                                                                  m_Delta_K_stack(t_Delta_K_stack),
                                                                                  m_Delta_dot_K_stack(t_Delta_dot_K_stack),
                                                                                  m_Delta_ddot_K_stack(t_Delta_ddot_K_stack),
                                                                                  m_angular_velocity_integrator(t_angular_velocity_integrator),
                                                                                  m_angular_acceleration_integrator(t_angular_acceleration_integrator),
                                                                                  m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {
        return m_Delta_ddot_K_stack->at(t_point)
                - ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ) * m_angular_acceleration_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_dot_K_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point);
    }

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_ddot_K_stack;


    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

};


struct DeltaLinearAccelerations : public OSNI::ODEAb {

    DeltaLinearAccelerations(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                             std::shared_ptr<const StrainParameterisationPerturbation> t_strain_parameterisation_perturbation,
                             std::shared_ptr<const IDMIntegrators> t_idm_integrators,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                             const unsigned int t_number_of_Chebyshev_points) :  OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                                                                                        t_number_of_Chebyshev_points),
                                                                                  m_K_stack(t_strain_parameterisation->m_K_stack),
                                                                                  m_dot_K_stack(t_strain_parameterisation->m_dot_K_stack),
                                                                                  m_Delta_K_stack(t_strain_parameterisation_perturbation->m_Delta_K_stack),
                                                                                  m_Delta_dot_K_stack(t_strain_parameterisation_perturbation->m_Delta_dot_K_stack),
                                                                                  m_Lambda_stack(t_strain_parameterisation->m_Lambda_stack),
                                                                                  m_dot_Lambda_stack(t_strain_parameterisation->m_dot_Lambda_stack),
                                                                                  m_Delta_Lambda_stack(t_strain_parameterisation_perturbation->m_Delta_Lambda_stack),
                                                                                  m_Delta_dot_Lambda_stack(t_strain_parameterisation_perturbation->m_Delta_dot_Lambda_stack),
                                                                                  m_Delta_ddot_Lambda_stack(t_strain_parameterisation_perturbation->m_Delta_ddot_Lambda_stack),
                                                                                  m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
                                                                                  m_angular_acceleration_integrator(t_idm_integrators->m_angular_acceleration),
                                                                                  m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
                                                                                  m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
                                                                                  m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity),
                                                                                  m_linear_acceleration_integrator(t_idm_integrators->m_linear_acceleration),
                                                                                  m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
        return -::LieAlgebra::skew( m_K_stack->at(t_point) );
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {
        return m_Delta_ddot_Lambda_stack->at(t_point)
                - ::LieAlgebra::skew( m_Lambda_stack->at(t_point) ) * m_Delta_angular_acceleration_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_Lambda_stack->at(t_point) ) * m_angular_acceleration_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ) * m_linear_acceleration_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_dot_Lambda_stack->at(t_point) ) * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_dot_K_stack->at(t_point) ) * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_dot_Lambda_stack->at(t_point) ) * m_angular_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_dot_K_stack->at(t_point) ) * m_linear_velocity_integrator->getStateAtPoint(t_point);
    }

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_dot_Lambda_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_ddot_Lambda_stack;


    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_acceleration_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;

};




struct DeltaInternalForcesIntegrator : public OSNI::ODEAb {

    DeltaInternalForcesIntegrator(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                                  std::shared_ptr<const StrainParameterisationPerturbation> t_strain_parameterisation_perturbation,
                                  std::shared_ptr<const IDMIntegrators> t_idm_integrators,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator,
                                  const Eigen::Matrix3d &t_Ml,
                                  const unsigned int t_number_of_Chebyshev_points) :    OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                                                                                                    t_number_of_Chebyshev_points),
                                                                                        m_K_stack(t_strain_parameterisation->m_K_stack),
                                                                                        m_Delta_K_stack(t_strain_parameterisation_perturbation->m_Delta_K_stack),
                                                                                        m_quaternion_integrator(t_idm_integrators->m_quaternion),
                                                                                        m_position_integrator(t_idm_integrators->m_position),
                                                                                        m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
                                                                                        m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
                                                                                        m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity),
                                                                                        m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
                                                                                        m_Delta_linear_acceleration_integrator(t_Delta_linear_acceleration_integrator),
                                                                                        m_Ml(t_Ml),
                                                                                        m_internal_forces_integrator(t_idm_integrators->m_internal_forces)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
        std::cout << "At point : " << t_point << "\n";
        std::cout << "  - A :\n" << -::LieAlgebra::skew( m_K_stack->at(t_point) ) << "\n";
        std::cout << " \n\n" ;

        return -::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {

        Eigen::Vector3d b = ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ).transpose()*m_internal_forces_integrator->getStateAtPoint(t_point)
                + m_Ml*m_Delta_linear_acceleration_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_Ml * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_Ml * m_linear_velocity_integrator->getStateAtPoint(t_point)
                - computeLocalExternalForces(m_quaternion_integrator->getStateAtPoint(t_point),
                                             m_position_integrator->getStateAtPoint(t_point));

        std::cout << "At point : " << t_point << "\n";
        std::cout << "  - b :\n" << b << "\n";
//        std::cout << "  - b :\n" << b << "\n";
//        std::cout << "  - b :\n" << b << "\n";
//        std::cout << "  - b :\n" << b << "\n";
//        std::cout << "  - b :\n" << b << "\n";
        std::cout << " \n\n" ;
        return b;

    }


    virtual Eigen::Vector3d computeLocalExternalForces(const Eigen::Vector4d &t_quaternion,
                                                       const Eigen::Vector3d &t_position) const
    {
        return Eigen::Vector3d::Zero();
    }





    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_position_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_acceleration_integrator;

    const Eigen::Matrix3d m_Ml;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;


};



struct TIDMIntegrators{

    TIDMIntegrators(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                    std::shared_ptr<const StrainParameterisationPerturbation> t_strain_parameterisation_perturbation,
                    std::shared_ptr<const IDMIntegrators> t_idm_integrators,
                    std::shared_ptr<const MaterialProperties> t_material_properties) :  m_strain_parameterisation(t_strain_parameterisation),
                                                                                        m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation),
                                                                                        m_idm_integrators(t_idm_integrators),
                                                                                        m_material_properties(t_material_properties)
    {}

    TIDMIntegrators(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                    std::shared_ptr<const StrainParameterisationPerturbation> t_strain_parameterisation_perturbation,
                    std::shared_ptr<const IDMIntegrators> t_idm_integrators,
                    std::shared_ptr<const MaterialProperties> t_material_properties,
                    const unsigned int t_number_of_Chebyshev_points) :  m_strain_parameterisation(t_strain_parameterisation),
                                                                        m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation),
                                                                        m_idm_integrators(t_idm_integrators),
                                                                        m_material_properties(t_material_properties),
                                                                        m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
    {}



    const std::shared_ptr<const StrainParameterisation> m_strain_parameterisation;
    const std::shared_ptr<const StrainParameterisationPerturbation> m_strain_parameterisation_perturbation;
    const std::shared_ptr<const IDMIntegrators> m_idm_integrators;

    const std::shared_ptr<const MaterialProperties> m_material_properties;

    const unsigned int m_number_of_Chebyshev_points { 17 };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_rotation { std::make_shared<DeltaRotation>( m_strain_parameterisation->m_K_stack,
                                                                                                  m_strain_parameterisation_perturbation->m_Delta_K_stack,
                                                                                                  m_number_of_Chebyshev_points) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_position { std::make_shared<DeltaPosition>( m_strain_parameterisation->m_K_stack,
                                                                                                  m_strain_parameterisation->m_Lambda_stack,
                                                                                                  m_strain_parameterisation_perturbation->m_Delta_Lambda_stack,
                                                                                                  m_Delta_rotation,
                                                                                                  m_number_of_Chebyshev_points) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_angular_velocity { std::make_shared<DeltaAngularVelocity>(m_strain_parameterisation->m_K_stack,
                                                                                                                m_strain_parameterisation_perturbation->m_Delta_K_stack,
                                                                                                                m_strain_parameterisation_perturbation->m_Delta_dot_K_stack,
                                                                                                                m_idm_integrators->m_angular_velocity,
                                                                                                                m_number_of_Chebyshev_points) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_linear_velocity { std::make_shared<DeltaLinearVelocity>(m_strain_parameterisation->m_K_stack,
                                                                                                              m_strain_parameterisation_perturbation->m_Delta_K_stack,
                                                                                                              m_strain_parameterisation->m_Lambda_stack,
                                                                                                              m_strain_parameterisation_perturbation->m_Delta_Lambda_stack,
                                                                                                              m_strain_parameterisation_perturbation->m_Delta_dot_Lambda_stack,
                                                                                                              m_idm_integrators->m_angular_velocity,
                                                                                                              m_Delta_angular_velocity,
                                                                                                              m_idm_integrators->m_linear_velocity,
                                                                                                              m_number_of_Chebyshev_points) };

    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_angular_acceleration { std::make_shared<DeltaAngularAccelerations>(m_strain_parameterisation->m_K_stack,
                                                                                                                         m_strain_parameterisation->m_dot_K_stack,
                                                                                                                         m_strain_parameterisation_perturbation->m_Delta_K_stack,
                                                                                                                         m_strain_parameterisation_perturbation->m_Delta_dot_K_stack,
                                                                                                                         m_strain_parameterisation_perturbation->m_Delta_ddot_K_stack,
                                                                                                                         m_idm_integrators->m_angular_velocity,
                                                                                                                         m_idm_integrators->m_angular_acceleration,
                                                                                                                         m_Delta_angular_velocity,
                                                                                                                         m_number_of_Chebyshev_points) };




    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_linear_acceleration { std::make_shared<DeltaLinearAccelerations>( m_strain_parameterisation,
                                                                                                                        m_strain_parameterisation_perturbation,
                                                                                                                        m_idm_integrators,
                                                                                                                        m_Delta_angular_velocity,
                                                                                                                        m_Delta_angular_acceleration,
                                                                                                                        m_Delta_linear_velocity,
                                                                                                                        m_number_of_Chebyshev_points) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_forces = { std::make_shared<DeltaInternalForcesIntegrator>(m_strain_parameterisation,
                                                                                                                          m_strain_parameterisation_perturbation,
                                                                                                                          m_idm_integrators,
                                                                                                                          m_Delta_angular_velocity,
                                                                                                                          m_Delta_linear_velocity,
                                                                                                                          m_Delta_linear_acceleration,
                                                                                                                          m_material_properties->m_M.block<3, 3>(3, 3),
                                                                                                                          m_number_of_Chebyshev_points) };

};

}   //  namespace CROSP

#endif // TIDM_INTEGRATORS_HPP
