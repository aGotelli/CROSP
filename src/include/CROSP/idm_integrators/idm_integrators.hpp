/**
 * \file idm_integrators.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the integrators for the Cosserat rod IDM
 * \date 13-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */



#ifndef IDM_INTEGRATORS_HPP
#define IDM_INTEGRATORS_HPP

#include <memory>

#include "CROSP/rod_properties/rod_properties.hpp"
#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include "OSNI/OSNI.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"


/// \brief CROSP::idm_integrators is the namespace containing the integrators for the IDM
namespace CROSP::idm_integrators {







struct ParameterisationStack {


    ParameterisationStack(const strain_parameterisation::StrainParameterisation &t_strain_parameterisation,
                          const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                          const unsigned int t_number_of_Chebyshev_points);


    /*!
     * \brief updateStacks updates the stacks of strain and its first and second derivative
     * \param t_qe the generalised elastic coordinates
     * \param t_dot_qe the first derivative of the generalised elastic coordinates
     * \param t_ddot_qe the second derivative of the generalised elastic coordinates
     */
    void updateStacks(const Eigen::VectorXd &t_qe,
                      const Eigen::VectorXd &t_dot_qe,
                      const Eigen::VectorXd &t_ddot_qe);



    //  The number of Chebyshev points used to represent the rod
    const unsigned int m_number_of_Chebyshev_points { 17 };


    //  The strain parameterisation of the rod
    strain_parameterisation::StrainParameterisation m_strain_parameterisation;

    //  The polynomial representation of the field of strain
    polynomial_representation::PolynomialRepresentation m_polynomial_representation;

    //  The stacks for the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Gamma_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };

    //  The stacks for the first derivative of the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Gamma_stack{
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };

    //  The stacks for the first derivative of the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Gamma_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };







    //  The vector stack of B Phi used to map the generalised coordinates into the strain field
    std::vector<Eigen::MatrixXd> m_map_to_strain_stack { [this](){

        std::vector<Eigen::MatrixXd> map_to_strain_stack(m_number_of_Chebyshev_points);

        std::vector<Eigen::MatrixXd> Phi_stack( m_number_of_Chebyshev_points );

        const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints( m_number_of_Chebyshev_points );

        std::generate(Phi_stack.begin(), Phi_stack.end(), [&, index=0]()mutable{
            const auto Phi = m_polynomial_representation.getPhi( Chebyshev_points[index] );
            index++;
            return Phi;
        });


        std::generate(map_to_strain_stack.begin(), map_to_strain_stack.end(), [&, index=0]()mutable{
            return m_polynomial_representation.m_B*Phi_stack[index++];
        });

        return map_to_strain_stack;
    }() };


};













/*!
 * \brief The QuaternionIntegrator struct implements the spectral numerical integration for the quaternion
 *
 * This class implemnts the integration for the quaternions forllowing the formula \f$ Q' = \frac{1}{2}  A_{(K)} Q
 */
struct QuaternionIntegrator : public OSNI::ODEA {

    QuaternionIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                         const double &t_upper_integration_limit=1.0f,
                         const Eigen::Vector4d &t_initial_condition=Eigen::Vector4d(1, 0, 0, 0));

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;

    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };

    Eigen::Matrix4d m_A_at_chebychev_point;

};



struct PositionIntegrator : public OSNI::ODEb {

    PositionIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                       std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                       const double &t_upper_integration_limit=1.0f,
                       const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_parameterisation_stack->m_Gamma_stack
    };

    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    Eigen::Quaterniond m_quaternion_at_point;

};




struct AngularVelocityIntegrator : public OSNI::ODEAb {

    AngularVelocityIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                              const double &t_upper_integration_limit=1.0f,
                              const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack {
        m_parameterisation_stack->m_dot_K_stack
    };


};


struct LinearVelocityIntegrator : public OSNI::ODEAb {

    LinearVelocityIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             const double &t_upper_integration_limit=1.0f,
                             const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point);


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point);


    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_parameterisation_stack->m_Gamma_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Gamma_stack {
        m_parameterisation_stack->m_dot_Gamma_stack
    };

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;



};




struct AngularAccelerationIntegrator : public OSNI::ODEAb {

    AngularAccelerationIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                  std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                  const double &t_upper_integration_limit=1.0f,
                                  const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack {
        m_parameterisation_stack->m_dot_K_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_K_stack {
        m_parameterisation_stack->m_ddot_K_stack
    };

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;


};







struct LinearAccelerationIntegrator : public OSNI::ODEAb {

    LinearAccelerationIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                                 std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                                 const double &t_upper_integration_limit=1.0f,
                                 const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_K_stack {
        m_parameterisation_stack->m_dot_K_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_parameterisation_stack->m_Gamma_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Gamma_stack{
        m_parameterisation_stack->m_dot_Gamma_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_Gamma_stack{
        m_parameterisation_stack->m_ddot_Gamma_stack
    };

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;



};




struct InternalForcesIntegrator : public OSNI::ODEAb {

    InternalForcesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                             std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_linear_acceleration_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                             std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                             const double &t_upper_integration_limit=1.0f,
                             const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());



    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computeDistributedForce(const unsigned int t_point) const;




    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };

    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    const Eigen::Matrix3d m_M_linear { m_rod_properties->getMLinear() };

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_acceleration;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_position;

};


struct InternalCouplesIntegrator : public OSNI::ODEAb {

    InternalCouplesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                              std::shared_ptr<const rod_properties::RodProperties> t_rod_properties,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_velocity_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_linear_velocity_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_angular_acceleration_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_quaternion_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_position_integrator,
                              std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                              const double &t_upper_integration_limit=1.0f,
                              const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero());




    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;


    virtual Eigen::VectorXd computeDistributedCouple(const unsigned int) const;


    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_K_stack {
        m_parameterisation_stack->m_K_stack
    };

    const std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_parameterisation_stack->m_Gamma_stack
    };

    const std::shared_ptr<const rod_properties::RodProperties> m_rod_properties;

    const Eigen::Matrix3d m_M_angular { m_rod_properties->getMAngular() };

    const Eigen::Matrix3d m_M_linear { m_rod_properties->getMLinear() };

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_linear_velocity;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_angular_acceleration;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_quaternion;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_position;

    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces;

};


struct GeneralisedForcesIntegrator : public OSNI::ODEb {

    GeneralisedForcesIntegrator(std::shared_ptr<const ParameterisationStack> t_parameterisation_stack,
                                std::shared_ptr<const OSNI::ODESolverInterface> t_internal_couples_integrator,
                                std::shared_ptr<const OSNI::ODESolverInterface> t_internal_forces_integrator,
                                const double &t_upper_integration_limit=1.0f);

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final;

    const std::shared_ptr<const ParameterisationStack> m_parameterisation_stack;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_couples_integrator;
    const std::shared_ptr<const OSNI::ODESolverInterface> m_internal_forces_integrator;

};









/*!
 * \brief The IDMIntegrators struct contains all the intgrators in order to integrate the kinematics and dynamics of the rod
 *
 * This struct is a collection of integrators to use in order to compute the Forward Kinematics and the Backward Dynamics
 */
struct IDMIntegrators {

    IDMIntegrators(const strain_parameterisation::StrainParameterisation t_strain_parameterisation,
                   const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                   const rod_properties::RodProperties t_rod_properties,
                   const unsigned int t_number_of_Chebyshev_points);


    //  Instance of the rod properties
    rod_properties::RodProperties m_rod_properties;

    std::shared_ptr<ParameterisationStack> m_parameterisation_stack;


//    //  Integrator for the quaternions
//    const std::shared_ptr<OSNI::ODESolverInterface> m_quaternion {
//        std::make_shared<QuaternionIntegrator>(m_parameterisation_stack,
//                                               m_rod_properties.m_rod_dimensions.m_L)
//    };

//    //  Integrator for the positions
//    const std::shared_ptr<OSNI::ODESolverInterface> m_position {
//        std::make_shared<PositionIntegrator>(m_parameterisation_stack,
//                                             m_quaternion,
//                                             m_rod_properties.m_rod_dimensions.m_L)
//    };

//    //  Integrator for the angular velocities
//    const std::shared_ptr<OSNI::ODESolverInterface> m_angular_velocity {
//        std::make_shared<AngularVelocityIntegrator>(m_parameterisation_stack,
//                                                    m_rod_properties.m_rod_dimensions.m_L)
//    };

//    //  Integrator for the linear velocities
//    const std::shared_ptr<OSNI::ODESolverInterface> m_linear_velocity {
//        std::make_shared<LinearVelocityIntegrator>(m_parameterisation_stack,
//                                                   m_angular_velocity,
//                                                   m_rod_properties.m_rod_dimensions.m_L)
//    };

//    //  Integrator for the angular accelerations
//    const std::shared_ptr<OSNI::ODESolverInterface> m_angular_acceleration {
//        std::make_shared<AngularAccelerationIntegrator>(m_parameterisation_stack,
//                                                        m_angular_velocity,
//                                                        m_rod_properties.m_rod_dimensions.m_L)
//    };

//    //  Integrator for the linear accelerations
//    const std::shared_ptr<OSNI::ODESolverInterface> m_linear_acceleration {
//        std::make_shared<LinearAccelerationIntegrator>(m_parameterisation_stack,
//                                                       m_angular_velocity,
//                                                       m_linear_velocity,
//                                                       m_angular_acceleration,
//                                                       m_rod_properties.m_rod_dimensions.m_L)
//    };


//    //  Integrator for the internal forces
//    const std::shared_ptr<OSNI::ODESolverInterface> m_internal_forces {
//        std::make_shared<InternalForcesIntegrator>(m_parameterisation_stack,
//                                                   m_rod_properties,
//                                                   m_angular_velocity,
//                                                   m_linear_velocity,
//                                                   m_linear_acceleration,
//                                                   m_quaternion,
//                                                   m_position,
//                                                   m_rod_properties.m_rod_dimensions.m_L)
//    };

//    //  Integrator for the internal couples
//    const std::shared_ptr<OSNI::ODESolverInterface> m_internal_couples {
//        std::make_shared<InternalCouplesIntegrator>(m_parameterisation_stack,
//                                                    m_rod_properties,
//                                                    m_angular_velocity,
//                                                    m_linear_velocity,
//                                                    m_angular_acceleration,
//                                                    m_quaternion,
//                                                    m_position,
//                                                    m_internal_forces,
//                                                    m_rod_properties.m_rod_dimensions.m_L)
//    };

//    //  Integrator for the generalised coordinates
//    const std::shared_ptr<OSNI::ODESolverInterface> m_generalised_forces {
//        std::make_unique<GeneralisedForcesIntegrator>(m_parameterisation_stack,
//                                                      m_internal_couples,
//                                                      m_internal_forces,
//                                                      m_rod_properties.m_rod_dimensions.m_L)
//    };


    /*!
     * \brief updateIntegrationDomain this function is used to set the integration upper limit at every integrator
     * \param t_upper_integration_limit is the upper integration limit used to scale the normalised domain [0, 1]
     */
    void updateIntegrationDomain(const double &t_upper_integration_limit=1.0f);
};


}   //  namespace CROSP::IDMIntegrators





#endif // IDM_INTEGRATORS_HPP
