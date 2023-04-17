#include "CROSP/CROSP/internal_actuation.hpp"


namespace CROSP::internal_actuation {




GeneralisedInternalActuationIntegrator::GeneralisedInternalActuationIntegrator(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation)
        : OSNI::ODEb(t_strain_parameterisation->m_polynomial_representation->getCoordinatesDimension(),
                     ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD,
                     t_strain_parameterisation->m_number_of_Chebyshev_points),
          m_strain_parameterisation(t_strain_parameterisation)
    {}

Eigen::VectorXd GeneralisedInternalActuationIntegrator::computerParametersVectorAtPoint(const unsigned int t_point)
{
    const double dc = 0.005;
    const Eigen::Vector3d D1(0,  dc, 0);
    const Eigen::Vector3d D2(0, -dc, 0);

    const auto D1_prime = Eigen::Vector3d::Zero();
    const auto D2_prime = Eigen::Vector3d::Zero();

    const Eigen::Vector3d Gamma = m_strain_parameterisation->m_Gamma_stack->at(t_point);
    const Eigen::Vector3d K = m_strain_parameterisation->m_K_stack->at(t_point);

    const auto K_skew = ::LieAlgebra::skew(K);

    const auto D1_skew = ::LieAlgebra::skew(D1);
    const auto D2_skew = ::LieAlgebra::skew(D2);

    const Eigen::Vector3d Gamma_1 = Gamma + K_skew * D1 + D1_prime;
    const Eigen::Vector3d Gamma_2 = Gamma + K_skew * D2 + D2_prime;

    const auto norm_Gamma_1 = Gamma_1.norm();
    const auto norm_Gamma_2 = Gamma_2.norm();

    ::LieAlgebra::Vector6d strain_map_1;
    strain_map_1 << D1_skew * Gamma_1   ,
                        Gamma_1         ;

    ::LieAlgebra::Vector6d strain_map_2;
    strain_map_2 << D2_skew * Gamma_2   ,
                        Gamma_2         ;


    const double tau_1 = m_tau(0);
    const double tau_2 = m_tau(1);

    const ::LieAlgebra::Vector6d strain_from_cable_1 = strain_map_1 * tau_1 / norm_Gamma_1;
    const ::LieAlgebra::Vector6d strain_from_cable_2 = strain_map_2 * tau_2 / norm_Gamma_2;

    const auto BPhi = m_strain_parameterisation->m_map_to_strain_stack.at(t_point);


    Eigen::VectorXd b = BPhi.transpose() * ( strain_from_cable_1 + strain_from_cable_2 );

    return b;
}




InternalActuation::InternalActuation(std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation)
    : m_integrator( GeneralisedInternalActuationIntegrator(t_strain_parameterisation) )
{}


void InternalActuation::setActuation([[maybe_unused]]const double &t_current_time)
{
    auto gain = Eigen::VectorXd(2);
    gain << 0.5,
            1;

    const double k = 30;
    const unsigned int number_of_cycles = 4;
    const double omega = 2*number_of_cycles*M_PI/(10);

    const auto theta = omega * t_current_time;

    m_integrator.m_tau = gain * k * (1 - cos( theta ) );
}

Eigen::VectorXd InternalActuation::computedGeneraliseInternalActuation()
{

    m_integrator.solveSystem();

    Eigen::VectorXd Q_ad = m_integrator.getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::BEGIN);

    return Q_ad;
}


}   //  namespace CROSP::internal_actuation
