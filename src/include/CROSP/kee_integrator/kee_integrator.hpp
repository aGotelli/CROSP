#ifndef KEE_INTEGRATOR_HPP
#define KEE_INTEGRATOR_HPP


#include <Eigen/Dense>

#include "OSNI/ODE/ode_b.hpp"




namespace CROSP::KeeIntegrator {



struct KeeIntegrator : public ::OSNI::ODEb {
    KeeIntegrator(const unsigned int t_coordinate_dimension,
                  const std::vector<Eigen::MatrixXd> &t_Phi_stack,
                  const Eigen::MatrixXd &t_Ha) :
        OSNI::ODEb(t_coordinate_dimension,
                   ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                   t_Phi_stack.size()),
        m_Phi_stack(t_Phi_stack),
        m_Ha(t_Ha)
    {}


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)final
    {
        const auto Phi = m_Phi_stack.at(t_point);

        return  Phi.transpose()*m_Ha*Phi;
    }

    const std::vector<Eigen::MatrixXd> m_Phi_stack;

    const Eigen::MatrixXd m_Ha;
};






}   //  namespace CROSP::KeeIntegrator

#endif // KEE_INTEGRATOR_HPP
