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



CosseratRod::CosseratRod(const strain_parameterisation::StrainParameterisation t_strain_parameterisation)
    : m_strain_parameterisation( t_strain_parameterisation )
{}



CosseratRod::CosseratRod(const strain_parameterisation::StrainParameterisation t_strain_parameterisation,
                         const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                         const rod_properties::RodDimensions t_rod_dimensions,
                         const rod_properties::MaterialProperties t_material_properties)
    : m_strain_parameterisation( t_strain_parameterisation ),
      m_polynomial_representation(t_polynomial_representation),
      m_rod_dimensions( t_rod_dimensions ),
      m_material_properties( t_material_properties )
{}




void CosseratRod::updateParameterisation(const Eigen::VectorXd &t_qe,
                                         const Eigen::VectorXd &t_dot_qe,
                                         const Eigen::VectorXd &t_ddot_qe)
{
    m_integrators->updateParameterisation(t_qe, t_dot_qe, t_ddot_qe);
}


void CosseratRod::forwardKinematics()
{
    m_integrators->forwardKinematics();
}


void CosseratRod::forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                    const Eigen::Vector3d &t_initial_position,
                                    const Eigen::Vector3d &t_initial_angular_velocity,
                                    const Eigen::Vector3d &t_initial_linear_velocity,
                                    const Eigen::Vector3d &t_initial_angular_acceleration,
                                    const Eigen::Vector3d &t_initial_linear_acceleration)
{
    m_integrators->forwardKinematics(t_initial_quaternion,
                                     t_initial_position,
                                     t_initial_angular_velocity,
                                     t_initial_linear_velocity,
                                     t_initial_angular_acceleration,
                                     t_initial_linear_acceleration);
}


void CosseratRod::forwardKinematics(const Eigen::Quaterniond &t_initial_quaternion,
                                    const Eigen::Vector3d &t_initial_position,
                                    const Eigen::Vector3d &t_initial_angular_velocity,
                                    const Eigen::Vector3d &t_initial_linear_velocity,
                                    const Eigen::Vector3d &t_initial_angular_acceleration,
                                    const Eigen::Vector3d &t_initial_linear_acceleration)
{
    const Eigen::Vector4d Q(t_initial_quaternion.w(),
                            t_initial_quaternion.x(),
                            t_initial_quaternion.y(),
                            t_initial_quaternion.z());

    forwardKinematics(Q, t_initial_position,
                      t_initial_angular_velocity, t_initial_linear_velocity,
                      t_initial_angular_acceleration, t_initial_linear_acceleration);
}






void CosseratRod::updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                              const Eigen::VectorXd &t_Delta_dot_qe,
                                              const Eigen::VectorXd &t_Delta_ddot_qe)
{

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
    m_integrators->updateDeltaParameterisation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);
}

void CosseratRod::forwardTangentKinematics()
{
    m_integrators->forwardTangentKinematics();
}


void CosseratRod::forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                           const Eigen::Vector3d &t_initial_Delta_position,
                                           const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                           const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                           const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                           const Eigen::Vector3d &t_initial_Delta_linear_acceleration)
{
    m_integrators->forwardTangentKinematics(t_initial_Delta_orientation,
                                            t_initial_Delta_position,
                                            t_initial_Delta_angular_velocity,
                                            t_initial_Delta_linear_velocity,
                                            t_initial_Delta_angular_acceleration,
                                            t_initial_Delta_linear_acceleration);
}






::LieAlgebra::Kinematics CosseratRod::getKinematicsAtTip()const
{
    return m_integrators->getKinematicsAtTip();
}


::LieAlgebra::TangentKinematics CosseratRod::getTangentKinematicsAtTip()const
{
    return m_integrators->getTangentKinematicsAtTip();
}



void CosseratRod::backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1)
{
    m_integrators->backwardDynamics(t_Lambda_X1);
}



void CosseratRod::backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
{
    m_integrators->backwardTangentDynamics(t_Delta_Lambda_X1);
}





::LieAlgebra::Vector6d CosseratRod::IDM(const Eigen::VectorXd &t_qe,
                                        const Eigen::VectorXd &t_dot_qe,
                                        const Eigen::VectorXd &t_ddot_qe,
                                        const Eigen::Vector4d &t_initial_quaternion,
                                        const Eigen::Vector3d &t_initial_position,
                                        const ::LieAlgebra::Vector6d &t_initial_twist,
                                        const ::LieAlgebra::Vector6d &t_initial_acceleration,
                                        const ::LieAlgebra::Vector6d &t_Lambda_X1)
{
    updateParameterisation(t_qe, t_dot_qe, t_ddot_qe);

    forwardKinematics(t_initial_quaternion,
                      t_initial_position,
                      t_initial_twist.block<3, 1>(0, 0),
                      t_initial_twist.block<3, 1>(3, 0),
                      t_initial_acceleration.block<3, 1>(0, 0),
                      t_initial_acceleration.block<3, 1>(3, 0));

    backwardDynamics(t_Lambda_X1);

    return getLambdaAtBase();
}


::LieAlgebra::Vector6d CosseratRod::IDM(const Eigen::VectorXd &t_qe,
                                        const Eigen::VectorXd &t_dot_qe,
                                        const Eigen::VectorXd &t_ddot_qe,
                                        const Eigen::Quaterniond &t_initial_quaternion,
                                        const Eigen::Vector3d &t_initial_position,
                                        const ::LieAlgebra::Vector6d &t_initial_twist,
                                        const ::LieAlgebra::Vector6d &t_initial_acceleration,
                                        const ::LieAlgebra::Vector6d &t_Lambda_X1)
{
    updateParameterisation(t_qe, t_dot_qe, t_ddot_qe);

    forwardKinematics(t_initial_quaternion,
                      t_initial_position,
                      t_initial_twist.block<3, 1>(0, 0),
                      t_initial_twist.block<3, 1>(3, 0),
                      t_initial_acceleration.block<3, 1>(0, 0),
                      t_initial_acceleration.block<3, 1>(3, 0));

    backwardDynamics(t_Lambda_X1);

    return getLambdaAtBase();
}

::LieAlgebra::Vector6d CosseratRod::TIDM(const Eigen::VectorXd &t_Delta_qe,
                                         const Eigen::VectorXd &t_Delta_dot_qe,
                                         const Eigen::VectorXd &t_Delta_ddot_qe,
                                         const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)
{
    updateDeltaParameterisation(t_Delta_qe, t_Delta_dot_qe, t_Delta_ddot_qe);


    forwardTangentKinematics();
    backwardTangentDynamics(t_Delta_Lambda_X1);

    return getDeltaLambdaAtBase();
}

::LieAlgebra::Vector6d CosseratRod::getLambdaAtBase()const
{
    return m_integrators->getLambdaAtBase();
}



LieAlgebra::Vector6d CosseratRod::getDeltaLambdaAtBase()const
{
    return m_integrators->getDeltaLambdaAtBase();
}




Eigen::VectorXd CosseratRod::getStaticInternalBalance(const Eigen::VectorXd &t_qe) const
{
    const Eigen::MatrixXd Kee = m_rod_properties.m_Kee;


    const Eigen::VectorXd Qe = Kee * t_qe;
    const Eigen::VectorXd Qa = m_integrators->getQaAtBase();
    const Eigen::VectorXd Q_ad = getQad();

    const Eigen::VectorXd internal_balance = Qe - Qa - Q_ad;
    return internal_balance;
}

Eigen::VectorXd CosseratRod::getInternalBalance(const Eigen::VectorXd &t_qe,
                                                const Eigen::VectorXd &t_dot_qe) const
{
    const Eigen::MatrixXd Dee = m_rod_properties.m_Dee;
    const Eigen::MatrixXd Kee = m_rod_properties.m_Kee;


    const Eigen::VectorXd Qe = Kee * t_qe;
    const Eigen::VectorXd Ce = Dee * t_dot_qe;
    const Eigen::VectorXd Qa = m_integrators->getQaAtBase();
    const Eigen::VectorXd Q_ad = getQad();

    const Eigen::VectorXd internal_balance = Qe + Ce - Qa - Q_ad;
    return internal_balance;

}





Eigen::VectorXd CosseratRod::getTangentStaticInternalBalance(const Eigen::VectorXd &t_Delta_qe)const
{
        Eigen::MatrixXd Kee = m_rod_properties.m_Kee;


    Eigen::VectorXd Delta_Qe = Kee * t_Delta_qe;
    Eigen::VectorXd Delta_Qa = m_integrators->getDeltaQaAtBase();

    Eigen::VectorXd Delta_internal_balance = Delta_Qe - Delta_Qa;

    return Delta_internal_balance;
}

Eigen::VectorXd CosseratRod::getTangentInternalBalance(const Eigen::VectorXd &t_Delta_qe,
                                                       const Eigen::VectorXd &t_Delta_dot_qe)const
{
    Eigen::MatrixXd Dee = m_rod_properties.m_Dee;
    Eigen::MatrixXd Kee = m_rod_properties.m_Kee;


    Eigen::VectorXd Delta_Qe = Kee * t_Delta_qe;
    Eigen::VectorXd Delta_Ce = Dee * t_Delta_dot_qe;
    Eigen::VectorXd Delta_Qa = m_integrators->getDeltaQaAtBase();

    Eigen::VectorXd Delta_internal_balance = Delta_Qe + Delta_Ce - Delta_Qa;
    return Delta_internal_balance;
}

void CosseratRod::updateInternalActuation([[maybe_unused]]const double &t_current_time) {};

Eigen::VectorXd CosseratRod::getQad()const
{
    return Eigen::VectorXd::Zero(m_polynomial_representation.getCoordinatesDimension());
}

Eigen::MatrixXd CosseratRod::getRodShapeFromElasticCoordinates(const Eigen::VectorXd &t_qe)const
{
    auto rod = *this;

    rod.updateParameterisation(t_qe, 0*t_qe, 0*t_qe);

    rod.forwardKinematics();

    return rod.getRodPositions();
}




void CosseratRod::printProperties()
{
    std::stringstream rod_properties;



    {   //  Start listing rod dimensions
    const auto l = m_rod_properties.m_rod_dimensions.m_L;
    rod_properties << "Rod dimensions :\n"
                      "     l : " << l << "\n" <<
                      m_rod_properties.m_rod_dimensions.m_cross_section->printProperties();
    }

    {   //  Start listing material properties
    const auto E = m_rod_properties.m_material_properties.m_E;
    const auto G = m_rod_properties.m_material_properties.m_G;
    const auto rho = m_rod_properties.m_material_properties.m_rho;
    rod_properties << "Material properties :\n"
                      "     E : " << std::setprecision(2) << std::scientific << E << "\n"
                      "     G : " << std::setprecision(2) << std::scientific << G << "\n"
                      "     ρ : " << std::setprecision(0) << std::fixed << rho << "\n";
    }



    {   //  Strain parameterisation
    rod_properties << "Strain parameterisation : \n";

    rod_properties << "     Number of Chebyshev points : " << m_integrators->m_number_of_Chebyshev_points << "\n";

       //  get the admitted deformations
    Eigen::VectorXi def(6);
    for(unsigned int i=0; i<6; i++)
        def[i] = m_polynomial_representation.m_admitted_deformations[i];
    rod_properties << "     Rod deformations : " << def.transpose() << "\n";


       //  Also print the number of modes
    Eigen::VectorXi ne_stack = Eigen::VectorXi::Zero(6);
    unsigned int j =0;
    for(unsigned int i=0; i<m_polynomial_representation.m_admitted_deformations.size(); i++)
        if(m_polynomial_representation.m_admitted_deformations[i])
            ne_stack[i] = m_polynomial_representation.m_number_of_modes_stack[j++];
    rod_properties << "     Number of modes  : " << ne_stack.transpose() << "\n";



       //  Details about the base
    const auto poly_base = m_polynomial_representation.m_polynomial_base;

    rod_properties << "     Polynomial base : ";

    std::string base = "Custom base";
    const unsigned int p = 5;
    const double x = 0.5;
    if( poly_base(p, x) == ::CROSP::polynomial_representation::legendre_polynomial_base(p, x) )
        base = "Legendre";

    if( poly_base(p, x) == ::CROSP::polynomial_representation::chebyshev_polynomial_base(p, x) )
        base = "Chebyshev";

    rod_properties << base;


    }

    std::cout << rod_properties.str() << "\n";
    std::cout.flush();
}

void CosseratRod::updateRodLength(const double &t_rod_lenght)
{
    m_rod_properties.updateRodLength( t_rod_lenght );

    m_integrators->updateIntegrationDomain( t_rod_lenght );
}

}   //  namespace CROSP
