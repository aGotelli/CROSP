/**
 * \file strain_parameterisation.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the functions and class for the strain parameterisation
 * \date 28-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */


#ifndef STRAIN_PARAMETERISATION_HPP
#define STRAIN_PARAMETERISATION_HPP

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <boost/math/special_functions/legendre.hpp>

namespace CROSP {



static Eigen::MatrixXd getPhi(const unsigned int t_ne,
                              const unsigned int t_na,
                              const double& t_X,
                              const std::function<double(const unsigned int, const double&)> t_polynomial_base=[](const unsigned int t_point, const double& t_x) {return boost::math::legendre_p(t_point, t_x);},
                              const double& t_begin=0,
                              const double& t_end=1)
{
    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    double x = ( 2 * t_X - ( t_end + t_begin) ) / ( t_end - t_begin );

    //  Compute the values of the polynomial for every element of the strain field
    Eigen::VectorXd Phi_i(t_ne, 1);
    for(unsigned int i=0; i<t_ne; i++)
        Phi_i[i] = t_polynomial_base(i, x);


    //  Define the matrix of bases
    Eigen::MatrixXd Phi = Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(t_na, t_na), Phi_i.transpose());


    return Phi;
}

static std::vector<Eigen::MatrixXd> generatePhiStack(const unsigned int t_ne,
                                                     const unsigned int t_na,
                                                     const std::vector<double> &t_Chebyshev_points,
                                                     const std::function<double(const unsigned int, const double&)> t_polynomial_base=[](const unsigned int t_point, const double& t_x) {return boost::math::legendre_p(t_point, t_x);})
{
    std::vector<Eigen::MatrixXd> Phi_stack( t_Chebyshev_points.size() );

    std::generate(Phi_stack.begin(), Phi_stack.end(), [&, index=0]()mutable{
        const auto Phi = getPhi(t_ne, t_na, t_Chebyshev_points[index], t_polynomial_base);
        index++;
        return Phi;
    });

    return Phi_stack;
}







struct StrainParameterisation {

    StrainParameterisation(const std::vector<Eigen::MatrixXd> &t_Phi_stack) : m_Phi_stack(t_Phi_stack) {}


    void update(const Eigen::VectorXd &t_qe,
                const Eigen::VectorXd &t_dot_qe,
                const Eigen::VectorXd &t_ddot_qe)
    {
        for(unsigned int i=0; i<m_stacks_dimension; i++){

            m_K_stack->at(i) = m_Phi_stack[i]*t_qe;
            m_dot_K_stack->at(i) = m_Phi_stack[i]*t_dot_qe;
            m_ddot_K_stack->at(i) = m_Phi_stack[i]*t_ddot_qe;

        }

        for(unsigned int i=0; i<m_stacks_dimension; i++){
            m_Lambda_stack->at(i) = Eigen::Vector3d(1, 0, 0);
            m_dot_Lambda_stack->at(i) = Eigen::Vector3d::Zero();
            m_ddot_Lambda_stack->at(i) = Eigen::Vector3d::Zero();
        }
    }

    const std::vector<Eigen::MatrixXd> m_Phi_stack;

    const unsigned int m_stacks_dimension { static_cast<unsigned int>( m_Phi_stack.size() ) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Lambda_stack{ std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };
};



}   //  namespace CROSP




#endif // STRAIN_PARAMETERISATION_HPP
