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

#include "CROSP/CROSP/cosserat_rod.hpp"


#include <boost/math/special_functions/legendre.hpp>

#include <string_view>
#include <iomanip>

namespace CROSP {

//template<class NumericalIntegrator>
//CosseratRod<NumericalIntegrator>::CosseratRod(polynomial_representation::PolynomialRepresentation t_polynomial_represenation,
//                                             unsigned int t_number_of_Chebyshev_points,
//                                             rod_properties::RodDimensions t_rod_dimensions,
//                                             rod_properties::MaterialProperties t_material_properties,
//                                             strain_parameterisation_stack::StrainFunction t_constrained_strain)
//    : m_polynomial_representation(t_polynomial_represenation),
//      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
//      m_constrained_strain(t_constrained_strain),
//      m_rod_properties(
//          std::make_shared<rod_properties::RodProperties>(t_rod_dimensions,
//                                                          t_material_properties)
//          )
//{}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::updateParameterisation(const Eigen::VectorXd &t_qe,
//                                         const Eigen::VectorXd &t_dot_qe,
//                                         const Eigen::VectorXd &t_ddot_qe)
//{
//    m_cosserat_rod_integrators->updateParameterisation(t_qe, t_dot_qe, t_ddot_qe);
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::forwardKinematics()
//{
//    m_cosserat_rod_integrators->forwardKinematics();
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
//                                    const Eigen::Vector3d &t_initial_position,
//                                    const Eigen::Vector3d &t_initial_angular_velocity,
//                                    const Eigen::Vector3d &t_initial_linear_velocity,
//                                    const Eigen::Vector3d &t_initial_angular_acceleration,
//                                    const Eigen::Vector3d &t_initial_linear_acceleration)
//{
//    m_cosserat_rod_integrators->forwardKinematics(t_initial_quaternion,
//                                                  t_initial_position,
//                                                  t_initial_angular_velocity,
//                                                  t_initial_linear_velocity,
//                                                  t_initial_angular_acceleration,
//                                                  t_initial_linear_acceleration);
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::forwardKinematics(const Eigen::Quaterniond &t_initial_quaternion,
//                                    const Eigen::Vector3d &t_initial_position,
//                                    const Eigen::Vector3d &t_initial_angular_velocity,
//                                    const Eigen::Vector3d &t_initial_linear_velocity,
//                                    const Eigen::Vector3d &t_initial_angular_acceleration,
//                                    const Eigen::Vector3d &t_initial_linear_acceleration)
//{
//    const Eigen::Vector4d Q(t_initial_quaternion.w(),
//                            t_initial_quaternion.x(),
//                            t_initial_quaternion.y(),
//                            t_initial_quaternion.z());

//    forwardKinematics(Q, t_initial_position,
//                      t_initial_angular_velocity, t_initial_linear_velocity,
//                      t_initial_angular_acceleration, t_initial_linear_acceleration);
//}






//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
//                                              const Eigen::VectorXd &t_Delta_dot_qe,
//                                              const Eigen::VectorXd &t_Delta_ddot_qe)
//{

//    if( std::find(t_Delta_qe.begin(), t_Delta_qe.end(), 1.0) == t_Delta_qe.end() or
//        std::accumulate(t_Delta_qe.begin(), t_Delta_qe.end(), 0) != 1.0){
//        std::stringstream help_message;
//        help_message << "You passed an invalid variation of the strain parameterisation. \n"
//                                "The function argument t_Delta_qe must contain only unitary entries and zeros. Here some examples : \n";
//        Eigen::VectorXd Delta_qe = Eigen::VectorXd::Zero(t_Delta_qe.size());
//        Delta_qe[0] = 1.0;
//        help_message << Delta_qe.transpose();
//        help_message << "\n" << "Or, again : \n";

//        Delta_qe.setZero();
//        Delta_qe[static_cast<unsigned int>(Delta_qe.size()/2.0)] = 1.0;
//        help_message << Delta_qe.transpose();

//        throw std::runtime_error{help_message.str()};

//    }


//    m_cosserat_rod_integrators->updateDeltaParameterisation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
//}

//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::forwardTangentKinematics()
//{
//    m_cosserat_rod_integrators->forwardTangentKinematics();
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
//                                           const Eigen::Vector3d &t_initial_Delta_position,
//                                           const Eigen::Vector3d &t_initial_Delta_angular_velocity,
//                                           const Eigen::Vector3d &t_initial_Delta_linear_velocity,
//                                           const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
//                                           const Eigen::Vector3d &t_initial_Delta_linear_acceleration)
//{
//    m_cosserat_rod_integrators->forwardTangentKinematics(t_initial_Delta_orientation,
//                                                         t_initial_Delta_position,
//                                                         t_initial_Delta_angular_velocity,
//                                                         t_initial_Delta_linear_velocity,
//                                                         t_initial_Delta_angular_acceleration,
//                                                         t_initial_Delta_linear_acceleration);
//}





//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//::LieAlgebra::Kinematics CosseratRod<NumericalIntegrator>::getKinematicsAtTip()const
//{
//    return m_cosserat_rod_integrators->getKinematicsAtTip();
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//::LieAlgebra::TangentKinematics CosseratRod<NumericalIntegrator>::getTangentKinematicsAtTip()const
//{
//    return m_cosserat_rod_integrators->getTangentKinematicsAtTip();
//}



//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1)
//{
//    m_cosserat_rod_integrators->backwardDynamics(t_Lambda_X1);
//}



//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
//{
//    m_cosserat_rod_integrators->backwardTangentDynamics(t_Delta_Lambda_X1);

//}




//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//::LieAlgebra::Vector6d CosseratRod<NumericalIntegrator>::IDM(const Eigen::VectorXd &t_qe,
//                                        const Eigen::VectorXd &t_dot_qe,
//                                        const Eigen::VectorXd &t_ddot_qe,
//                                        const Eigen::Vector4d &t_initial_quaternion,
//                                        const Eigen::Vector3d &t_initial_position,
//                                        const ::LieAlgebra::Vector6d &t_initial_twist,
//                                        const ::LieAlgebra::Vector6d &t_initial_acceleration,
//                                        const ::LieAlgebra::Vector6d &t_Lambda_X1)
//{
//    m_cosserat_rod_integrators->updateParameterisation(t_qe, t_dot_qe, t_ddot_qe);

//    m_cosserat_rod_integrators->forwardKinematics(t_initial_quaternion,
//                                                  t_initial_position,
//                                                  t_initial_twist.block<3, 1>(0, 0),
//                                                  t_initial_twist.block<3, 1>(3, 0),
//                                                  t_initial_acceleration.block<3, 1>(0, 0),
//                                                  t_initial_acceleration.block<3, 1>(3, 0));

//    m_cosserat_rod_integrators->backwardDynamics(t_Lambda_X1);

//    return getLambdaAtBase();
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//::LieAlgebra::Vector6d CosseratRod<NumericalIntegrator>::IDM(const Eigen::VectorXd &t_qe,
//                                        const Eigen::VectorXd &t_dot_qe,
//                                        const Eigen::VectorXd &t_ddot_qe,
//                                        const Eigen::Quaterniond &t_initial_quaternion,
//                                        const Eigen::Vector3d &t_initial_position,
//                                        const ::LieAlgebra::Vector6d &t_initial_twist,
//                                        const ::LieAlgebra::Vector6d &t_initial_acceleration,
//                                        const ::LieAlgebra::Vector6d &t_Lambda_X1)
//{
//    const Eigen::Vector4d initial_quaternion(t_initial_quaternion.w(),
//                                             t_initial_quaternion.x(),
//                                             t_initial_quaternion.y(),
//                                             t_initial_quaternion.z());

//    return IDM(t_qe,
//               t_dot_qe,
//               t_ddot_qe,
//               initial_quaternion,
//               t_initial_position,
//               t_initial_twist,
//               t_initial_acceleration,
//               t_Lambda_X1);
//}

//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//::LieAlgebra::Vector6d CosseratRod<NumericalIntegrator>::TIDM(const Eigen::VectorXd &t_Delta_qe,
//                                         const Eigen::VectorXd &t_Delta_dot_qe,
//                                         const Eigen::VectorXd &t_Delta_ddot_qe,
//                                         const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
//{
//    updateDeltaParameterisation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);


//    forwardTangentKinematics();
//    backwardTangentDynamics(t_Delta_Lambda_X1);

//    return getDeltaLambdaAtBase();
//}

//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//::LieAlgebra::Vector6d CosseratRod<NumericalIntegrator>::getLambdaAtBase()const
//{
//    return m_cosserat_rod_integrators->getLambdaAtBase();
//}



//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//::LieAlgebra::Vector6d CosseratRod<NumericalIntegrator>::getDeltaLambdaAtBase()const
//{
//    return m_cosserat_rod_integrators->getDeltaLambdaAtBase();
//}



//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//Eigen::VectorXd CosseratRod<NumericalIntegrator>::getStaticInternalBalance(const Eigen::VectorXd &t_qe) const
//{

//    const Eigen::VectorXd Qe = m_Kee * t_qe;
//    const Eigen::VectorXd Qa = m_cosserat_rod_integrators->getQaAtBase();
//    const Eigen::VectorXd Q_ad = m_cosserat_rod_integrators->getQad();

//    const Eigen::VectorXd internal_balance = Qe - Qa - Q_ad;
//    return internal_balance;
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//Eigen::VectorXd CosseratRod<NumericalIntegrator>::getInternalBalance(const Eigen::VectorXd &t_qe,
//                                                const Eigen::VectorXd &t_dot_qe) const
//{

//    const Eigen::VectorXd Qe = m_Dee * t_qe;
//    const Eigen::VectorXd Ce = m_Kee * t_dot_qe;
//    const Eigen::VectorXd Qa = m_cosserat_rod_integrators->getQaAtBase();
//    const Eigen::VectorXd Q_ad = m_cosserat_rod_integrators->getQad();

//    const Eigen::VectorXd internal_balance = Qe + Ce - Qa - Q_ad;
//    return internal_balance;

//}





//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//Eigen::VectorXd CosseratRod<NumericalIntegrator>::getTangentStaticInternalBalance(const Eigen::VectorXd &t_Delta_qe)const
//{

//    Eigen::VectorXd Delta_Qe = m_Kee * t_Delta_qe;
//    Eigen::VectorXd Delta_Qa = m_cosserat_rod_integrators->getDeltaQaAtBase();

//    Eigen::VectorXd Delta_internal_balance = Delta_Qe - Delta_Qa;

//    return Delta_internal_balance;
//}



//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//Eigen::VectorXd CosseratRod<NumericalIntegrator>::getTangentInternalBalance(const Eigen::VectorXd &t_Delta_qe,
//                                                       const Eigen::VectorXd &t_Delta_dot_qe)const
//{

//    Eigen::VectorXd Delta_Qe = m_Kee * t_Delta_qe;
//    Eigen::VectorXd Delta_Ce = m_Dee * t_Delta_dot_qe;
//    Eigen::VectorXd Delta_Qa = m_cosserat_rod_integrators->getDeltaQaAtBase();


//    Eigen::VectorXd Delta_internal_balance = Delta_Qe + Delta_Ce - Delta_Qa;
//    return Delta_internal_balance;
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::updateInternalActuation(const double &t_current_time)
//{
//    m_cosserat_rod_integrators->updateInternalActuation(t_current_time);
//}



//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//Eigen::MatrixXd CosseratRod<NumericalIntegrator>::getRodShapeFromElasticCoordinates(const Eigen::VectorXd &t_qe)
//{

//    updateParameterisation(t_qe, 0*t_qe, 0*t_qe);

//    forwardKinematics();

//    return getRodPositionsAtChebyshevPoints();
//}




//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::printProperties()
//{
//    std::stringstream rod_properties;



//    {   //  Start listing rod dimensions
//    const auto l = m_rod_properties->m_rod_dimensions.m_L;
//    rod_properties << "Rod dimensions :\n"
//                      "     l : " << l << "\n" <<
//                      m_rod_properties->m_rod_dimensions.m_cross_section->printProperties();
//    }

//    {   //  Start listing material properties
//    const auto E = m_rod_properties->m_material_properties.m_E;
//    const auto G = m_rod_properties->m_material_properties.m_G;
//    const auto rho = m_rod_properties->m_material_properties.m_rho;
//    rod_properties << "Material properties :\n"
//                      "     E : " << std::setprecision(2) << std::scientific << E << "\n"
//                      "     G : " << std::setprecision(2) << std::scientific << G << "\n"
//                      "     ρ : " << std::setprecision(0) << std::fixed << rho << "\n";
//    }



//    {   //  Strain parameterisation
//    rod_properties << "Strain parameterisation : \n";


//       //  get the admitted deformations
//    Eigen::VectorXi def(6);
//    for(unsigned int i=0; i<6; i++)
//        def[i] = m_polynomial_representation.m_admitted_deformations[i];
//    rod_properties << "     Rod deformations : " << def.transpose() << "\n";


//       //  Also print the number of modes
//    Eigen::VectorXi ne_stack = Eigen::VectorXi::Zero(6);
//    unsigned int j =0;
//    for(unsigned int i=0; i<m_polynomial_representation.m_admitted_deformations.size(); i++)
//        if(m_polynomial_representation.m_admitted_deformations[i])
//            ne_stack[i] = m_polynomial_representation.m_number_of_modes_stack[j++];
//    rod_properties << "     Number of modes  : " << ne_stack.transpose() << "\n";



//       //  Details about the base
//    const auto poly_base = m_polynomial_representation.m_polynomial_base;

//    rod_properties << "     Polynomial base : ";

//    std::string base = "Custom base";
//    const unsigned int p = 5;
//    const double x = 0.5;
//    if( poly_base(p, x) == ::CROSP::polynomial_representation::legendre_polynomial_base(p, x) )
//        base = "Legendre";

//    if( poly_base(p, x) == ::CROSP::polynomial_representation::chebyshev_polynomial_base(p, x) )
//        base = "Chebyshev";

//    rod_properties << base;


//    }

//    std::cout << rod_properties.str() << "\n";
//    std::cout.flush();
//}

//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::updateRodLength(const double &t_rod_lenght)
//{
//    const double scale = t_rod_lenght/m_rod_properties->m_rod_dimensions.m_L;
//    m_rod_properties->m_rod_dimensions.m_L = t_rod_lenght;

//    m_Kee *= scale;
//    m_Dee *= scale;

//    m_cosserat_rod_integrators->updateIntegrationDomain(t_rod_lenght);
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::updateRodProperties(const rod_properties::RodDimensions &t_rod_dimensions,
//                                           const rod_properties::MaterialProperties &t_material_properties)
//{
////    m_rod_properties->updateRodProperties(t_rod_dimensions,
////                                          t_material_properties,
////                                          m_strain_parameterisation->m_polynomial_representation);
//}


//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//void CosseratRod<NumericalIntegrator>::updateRodProperties(const double &t_EI)
//{
//    m_rod_properties->m_H(2, 2) = t_EI;

////    m_Kee = defineKee();

//    m_Dee = m_rod_properties->m_material_properties.m_mu*m_Kee;
//}



//template<numerical_integrators::CosseratIntegrator NumericalIntegrator>
//Eigen::MatrixXd CosseratRod<NumericalIntegrator>::defineKee()
//{

//    const unsigned int n = m_polynomial_representation.getCoordinatesDimension();

//    const Eigen::MatrixXd Ha = m_polynomial_representation.m_B.transpose() * m_rod_properties->m_H * m_polynomial_representation.m_B;


//    typedef boost::numeric::odeint::runge_kutta_dopri5< Eigen::MatrixXd, double,
//                                                         Eigen::MatrixXd, double,
//                                                         boost::numeric::odeint::vector_space_algebra> Ke_stepper;
//    Eigen::MatrixXd Kee = Eigen::MatrixXd::Zero(n, n);

//    const double X0 = 0.0;
//    const double X1 = 1.0;
//    const double dX = 0.0005;

//    boost::numeric::odeint::integrate_adaptive(Ke_stepper(), [&](const Eigen::MatrixXd &, Eigen::MatrixXd &t_dKeeds, const double t_X){
//        const auto Phi = m_polynomial_representation.getPhi( t_X );

//        t_dKeeds = Phi.transpose()*Ha*Phi;
//    }, Kee, X0, X1, dX);

//    Kee *= m_rod_properties->m_rod_dimensions.m_L;

//    return Kee;
//}


}   //  namespace CROSP
