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
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_acceleration_integrator,
                                  std::shared_ptr<const MaterialProperties> t_material_properties,
                                  const unsigned int t_number_of_Chebyshev_points) :    OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                                                                                                    t_number_of_Chebyshev_points),
                                                                                        m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
                                                                                        m_K_stack(t_strain_parameterisation->m_K_stack),
                                                                                        m_Delta_K_stack(t_strain_parameterisation_perturbation->m_Delta_K_stack),
                                                                                        m_quaternion_integrator(t_idm_integrators->m_quaternion),
                                                                                        m_position_integrator(t_idm_integrators->m_position),
                                                                                        m_Delta_rotation_integrator(t_Delta_rotation_integrator),
                                                                                        m_angular_velocity_integrator(t_idm_integrators->m_angular_velocity),
                                                                                        m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
                                                                                        m_linear_velocity_integrator(t_idm_integrators->m_linear_velocity),
                                                                                        m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
                                                                                        m_Delta_linear_acceleration_integrator(t_Delta_linear_acceleration_integrator),
                                                                                        m_Ml(t_material_properties->m_M.block<3,3>(3, 3)),
                                                                                        m_material_properties(t_material_properties),
                                                                                        m_internal_forces_integrator(t_idm_integrators->m_internal_forces)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
//        std::cout << "At point : " << t_point << "\n";
//        std::cout << "  - A :\n" << -::LieAlgebra::skew( m_K_stack->at(t_point) ) << "\n";
//        std::cout << " \n\n" ;

        return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {

//        std::cout << "At point : " << t_point << "\n";
//        std::cout << "  - Omega : \n" << m_angular_velocity_integrator->getStateAtPoint(t_point) << "\n";
//        std::cout << "  - V : \n" << m_linear_velocity_integrator->getStateAtPoint(t_point) << "\n";
//        std::cout << "  - Delta Omega : \n" << m_Delta_angular_velocity_integrator->getStateAtPoint(t_point) << "\n";
//        std::cout << "  - Delta V : \n" << m_Delta_linear_velocity_integrator->getStateAtPoint(t_point) << "\n";
//        std::cout << "  - Delta dot V : \n" << m_Delta_linear_acceleration_integrator->getStateAtPoint(t_point) << "\n";
//        std::cout << "  - N : \n" << m_internal_forces_integrator->getStateAtPoint(t_point) << "\n";



        const Eigen::Vector3d Delta_N_ext = computeLocalExternalForces(t_point);
        const Eigen::Vector3d Delta_K_Lambda = ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ).transpose()*m_internal_forces_integrator->getStateAtPoint(t_point);
        const Eigen::Vector3d M_Delta_dot_eta = m_Ml*m_Delta_linear_acceleration_integrator->getStateAtPoint(t_point);
        const Eigen::Vector3d m_Omega_M_Delta_V = ::LieAlgebra::skew( m_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_Ml * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point);
        const Eigen::Vector3d m_Delta_Omega_M_V = ::LieAlgebra::skew( m_Delta_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_Ml * m_linear_velocity_integrator->getStateAtPoint(t_point);

//        std::cout << "  - M_Delta_dot_eta : \n" << M_Delta_dot_eta << "\n";
//        std::cout << "  - m_Delta_Omega_M_V : \n" << m_Delta_Omega_M_V << "\n";
//        std::cout << "  - m_Omega_M_Delta_V : \n" << m_Omega_M_Delta_V << "\n";
//        std::cout << "  - Delta_K_Lambda : \n" << Delta_K_Lambda << "\n";
//        std::cout << "  - Delta_N_ext : \n" << Delta_N_ext << "\n";




        Eigen::Vector3d b = Eigen::Vector3d::Zero()
                + m_Ml*m_Delta_linear_acceleration_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_Delta_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_Ml * m_linear_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew( m_angular_velocity_integrator->getStateAtPoint(t_point) ).transpose() * m_Ml * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
                - computeLocalExternalForces(t_point)
                + ::LieAlgebra::skew( m_Delta_K_stack->at(t_point) ).transpose()*m_internal_forces_integrator->getStateAtPoint(t_point)
                ;


//        std::cout << "  - b :\n" << b << "\n";

//        std::cout << "\n\n\n";

        b_stack(Eigen::all, t_point) = b;

//        std::cout << b_stack << std::endl;
//                std::cout << "\n\n\n";


        return b;

    }


    virtual Eigen::Vector3d computeLocalExternalForces(unsigned int t_point) const
    {
        if(t_point == 0 or t_point == m_number_of_Chebyshev_points-1)
            return Eigen::Vector3d::Zero();

        Eigen::Vector4d q(m_quaternion_integrator->getStateAtPoint(t_point));

        Eigen::Matrix3d R = Eigen::Quaterniond(q(0),
                                               q(1),
                                               q(2),
                                               q(3)).toRotationMatrix();

        Eigen::Vector3d Delta_rotation = m_Delta_rotation_integrator->getStateAtPoint(t_point);

        Eigen::Vector3d Delta_N_bar = -::LieAlgebra::skew( Delta_rotation ).transpose()
                                        *R.transpose()*m_material_properties->m_rho*m_material_properties->m_A*Eigen::Vector3d(0, 0, 9.81);


        return Delta_N_bar;
    }

//    void printPointInfo(unsigned int t_point, unsigned int t_number_of_Chebyshev_points)
//    {
//        std::cout << "At point : " << t_point << ", corresponding to X = " << ::Chebyshev::ComputeChebyshevPoints(t_number_of_Chebyshev_points)[t_point] << "\n";
//        std::cout
//    }

    const unsigned int m_number_of_Chebyshev_points;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack;
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_position_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_rotation_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_acceleration_integrator;

    const Eigen::Matrix3d m_Ml;
    const std::shared_ptr<const MaterialProperties> m_material_properties;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;

    Eigen::MatrixXd b_stack { Eigen::MatrixXd::Ones(3, /*m_K_stack->size()*/m_quaternion_integrator->getStackAsMatrix().cols()) };


};






struct DeltaInternalCouplesIntegrator : public OSNI::ODEAb {

    DeltaInternalCouplesIntegrator(std::shared_ptr<const MaterialProperties> t_material_properties,
                                   std::shared_ptr<const IDMIntegrators> t_idm_integrators,
                                   std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                                   std::shared_ptr<const StrainParameterisationPerturbation> t_strain_parameterisation_perturbation,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_rotation_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_angular_acceleration_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_linear_velocity_integrator,
                                   std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                   const unsigned int t_number_of_Chebyshev_points) :    OSNI::ODEAb(3, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                                                                                                    t_number_of_Chebyshev_points),
                                                                                        m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
                                                                                        m_material_properties(t_material_properties),
                                                                                        m_idm_integrators(t_idm_integrators),
                                                                                        m_strain_parameterisation(t_strain_parameterisation),
                                                                                        m_strain_parameterisation_perturbation(t_strain_parameterisation_perturbation),
                                                                                        m_Delta_rotation_integrator(t_Delta_rotation_integrator),
                                                                                        m_Delta_angular_velocity_integrator(t_Delta_angular_velocity_integrator),
                                                                                        m_Delta_angular_acceleration_integrator(t_Delta_angular_acceleration_integrator),
                                                                                        m_Delta_linear_velocity_integrator(t_Delta_linear_velocity_integrator),
                                                                                        m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator)
    {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)
    {
        return ::LieAlgebra::skew( m_K_stack->at(t_point) ).transpose();
    }

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)
    {

        Eigen::Vector3d b = Eigen::Vector3d::Zero()
                + ::LieAlgebra::skew(m_Lambda_stack->at(t_point)).transpose() * m_Delta_internal_forces_integrator->getStateAtPoint(t_point)
                + ::LieAlgebra::skew(m_Delta_K_stack->at(t_point)).transpose() * m_internal_couples_integrator->getStateAtPoint(t_point)
                + ::LieAlgebra::skew(m_Delta_Lambda_stack->at(t_point)).transpose() * m_internal_forces_integrator->getStateAtPoint(t_point)
                + m_M_angular * m_Delta_angular_acceleration_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew(m_angular_velocity_integrator->getStateAtPoint(t_point)).transpose() * m_M_angular * m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew(m_linear_velocity_integrator->getStateAtPoint(t_point)).transpose() * m_M_linear * m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew(m_Delta_angular_velocity_integrator->getStateAtPoint(t_point)).transpose() * m_M_angular * m_angular_velocity_integrator->getStateAtPoint(t_point)
                - ::LieAlgebra::skew(m_Delta_linear_velocity_integrator->getStateAtPoint(t_point)).transpose() * m_M_linear * m_linear_velocity_integrator->getStateAtPoint(t_point)
                - computeLocalExternalCouples(t_point)
                ;

        b_stack(Eigen::all, t_point) = b;

        return b;

    }


    virtual Eigen::Vector3d computeLocalExternalCouples(unsigned int t_point) const
    {
        return Eigen::Vector3d::Zero();
    }



    const unsigned int m_number_of_Chebyshev_points;

    const std::shared_ptr<const MaterialProperties> m_material_properties;

    const std::shared_ptr<const IDMIntegrators> m_idm_integrators;

    const std::shared_ptr<const StrainParameterisation> m_strain_parameterisation;

    const std::shared_ptr<const StrainParameterisationPerturbation> m_strain_parameterisation_perturbation;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion_integrator { m_idm_integrators->m_quaternion };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_position_integrator { m_idm_integrators->m_position };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack { m_strain_parameterisation->m_K_stack };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_K_stack { m_strain_parameterisation_perturbation->m_Delta_K_stack };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Lambda_stack { m_strain_parameterisation->m_Lambda_stack };
    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack { m_strain_parameterisation_perturbation->m_Delta_Lambda_stack};



    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_rotation_integrator;


    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity_integrator { m_idm_integrators->m_angular_velocity };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_velocity_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_angular_acceleration_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity_integrator { m_idm_integrators->m_linear_velocity };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_linear_velocity_integrator;


    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator { m_idm_integrators->m_internal_forces };
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_forces_integrator;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator { m_idm_integrators->m_internal_couples };


    const Eigen::Matrix3d m_M_angular { m_material_properties->m_M.block<3,3>(0, 0) };
    const Eigen::Matrix3d m_M_linear { m_material_properties->m_M.block<3,3>(3, 3) };


    Eigen::MatrixXd b_stack { Eigen::MatrixXd::Ones(3, /*m_K_stack->size()*/m_quaternion_integrator->getStackAsMatrix().cols()) };


};



struct DeltaGeneralisedForcesIntegrator : public OSNI::ODEb {

    DeltaGeneralisedForcesIntegrator(std::shared_ptr<const StrainParameterisation> t_strain_parameterisation,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_couples_integrator,
                                     std::shared_ptr<const OSNI::ODESolverInterface> t_Delta_internal_forces_integrator,
                                     const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEb(t_strain_parameterisation->m_ne,
                                                                                                   ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                                                                                                   t_number_of_Chebyshev_points),
                                                                                        m_strain_parameterisation(t_strain_parameterisation),
                                                                                        m_Delta_internal_couples_integrator(t_Delta_internal_couples_integrator),
                                                                                        m_Delta_internal_forces_integrator(t_Delta_internal_forces_integrator)
    {}

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        Eigen::Vector3d Delta_C = m_Delta_internal_couples_integrator->getStateAtPoint(t_point);
        Eigen::Vector3d Delta_N = m_Delta_internal_forces_integrator->getStateAtPoint(t_point);

        Eigen::Matrix<double, 6, 1> Delta_Lambda;
        Delta_Lambda << Delta_C,
                        Delta_N;

        const auto Phi = m_strain_parameterisation->m_Phi_stack.at(t_point);
        const auto B = m_strain_parameterisation->m_B;
        return -Phi.transpose()*B.transpose()*Delta_Lambda;
    }


    const std::shared_ptr<const StrainParameterisation> m_strain_parameterisation;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_couples_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_Delta_internal_forces_integrator;

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


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_forces { std::make_shared<DeltaInternalForcesIntegrator>(m_strain_parameterisation,
                                                                                                                        m_strain_parameterisation_perturbation,
                                                                                                                        m_idm_integrators,
                                                                                                                        m_Delta_rotation,
                                                                                                                        m_Delta_angular_velocity,
                                                                                                                        m_Delta_linear_velocity,
                                                                                                                        m_Delta_linear_acceleration,
                                                                                                                        m_material_properties,
                                                                                                                        m_number_of_Chebyshev_points) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_internal_couples { std::make_shared<DeltaInternalCouplesIntegrator>(m_material_properties,
                                                                                                                          m_idm_integrators,
                                                                                                                          m_strain_parameterisation,
                                                                                                                          m_strain_parameterisation_perturbation,
                                                                                                                          m_Delta_rotation,
                                                                                                                          m_Delta_angular_velocity,
                                                                                                                          m_Delta_angular_acceleration,
                                                                                                                          m_Delta_linear_velocity,
                                                                                                                          m_Delta_internal_forces,
                                                                                                                          m_number_of_Chebyshev_points) };


    std::shared_ptr<OSNI::ODESolverInterface> m_Delta_generalised_forces { std::make_shared<DeltaGeneralisedForcesIntegrator>(m_strain_parameterisation,
                                                                                                                              m_Delta_internal_couples,
                                                                                                                              m_Delta_internal_forces,
                                                                                                                              m_number_of_Chebyshev_points) };





};

}   //  namespace CROSP

#endif // TIDM_INTEGRATORS_HPP
