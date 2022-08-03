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
#include <algorithm>

#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <boost/math/special_functions/legendre.hpp>

#include "material_properties.hpp"

#include "math_tools/Chebyshev/chebyshev_differentiation.hpp"
#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

namespace CROSP {

namespace StrainParameterisation {


typedef std::function<double(const unsigned int, const double&)> BaseFunction;
static const BaseFunction default_base_function { [](const unsigned int t_point,
                                                     const double& t_x)
                                                        {
                                                            return boost::math::legendre_p(t_point, t_x);
                                                        }
                                                };


Eigen::MatrixXd getPhi(const unsigned int t_ne,
                              const unsigned int t_na,
                              const double& t_X,
                              const BaseFunction t_polynomial_base=default_base_function,
                              const double& t_begin=0,
                              const double& t_end=1);

std::vector<Eigen::MatrixXd> generatePhiStack(const unsigned int t_ne,
                                              const unsigned int t_na,
                                              const std::vector<double> &t_Chebyshev_points,
                                              const BaseFunction t_polynomial_base=default_base_function);







struct StrainParameterisation {

    StrainParameterisation()=default;


    StrainParameterisation(unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const unsigned int t_number_of_Chebyshev_points,
                           const Eigen::VectorXd &t_constrained_strain,
                           const BaseFunction t_polynomial_base=default_base_function);


    inline unsigned int getNumberOfChebyshewPoints()const{return m_number_of_Chebyshev_points;}


    void update(const Eigen::VectorXd &t_qe,
                const Eigen::VectorXd &t_dot_qe,
                const Eigen::VectorXd &t_ddot_qe);




    const std::vector<bool> m_admitted_deformations { true, true, true, false, false, false };

    const unsigned int m_na { [&]()->unsigned int{ return std::count(m_admitted_deformations.begin(), m_admitted_deformations.end(), true);}() };

    const unsigned int m_ne { 3 };

    const Eigen::MatrixXd m_B { [&](){
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);

            std::vector<int> indexes;
            std::for_each(m_admitted_deformations.begin(),
                          m_admitted_deformations.end(),
                          [&indexes, index=0](const bool dof)mutable{   if(dof == true)
                                                                            indexes.push_back(index);
                                                                        index++;});
            const Eigen::MatrixXd map = I(Eigen::all, indexes);
            return map;
        }() };

    const Eigen::MatrixXd m_B_bar { [&](){

            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);

            std::vector<int> indexes;
            std::for_each(m_admitted_deformations.begin(),
                          m_admitted_deformations.end(),
                          [&indexes, index=0](const bool dof)mutable{   if(dof == false)
                                                                            indexes.push_back(index);
                                                                        index++;});
            const Eigen::MatrixXd map = I(Eigen::all, indexes);
            return map;
        }() };

    const Eigen::VectorXd m_constrained_strain { [&](){

            //  constrain only the evolution along x
            Eigen::VectorXd constrains(6);
            constrains <<   0,
                            0,
                            0,
                            1.0,
                            0,
                            0;

            std::vector<int> indexes;
            std::for_each(m_admitted_deformations.begin(),
                          m_admitted_deformations.end(),
                          [&indexes, index=0](const bool dof)mutable{   if(dof == false)
                                                                            indexes.push_back(index);
                                                                        index++;});
            Eigen::VectorXd xi_c = constrains(indexes);

            return xi_c;
        }() };



    const unsigned int m_number_of_Chebyshev_points { 17 };

    const std::function<double(const unsigned int, const double&)> m_polynomial_base { [](const unsigned int t_point, const double& t_x) {return boost::math::legendre_p(t_point, t_x);} };

    const std::vector<Eigen::MatrixXd> m_Phi_stack { generatePhiStack(m_ne, m_na, ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points), m_polynomial_base) };

    const std::vector<Eigen::MatrixXd> m_strains_map_stack { [&](){
            std::vector<Eigen::MatrixXd> strains_map_stack(m_number_of_Chebyshev_points);

            std::generate(strains_map_stack.begin(), strains_map_stack.end(), [&, index=0]()mutable{
               return m_B*m_Phi_stack[index++];
            });

            return strains_map_stack;
        }() };



    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Lambda_stack{ std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };
};


struct StrainParameterisationPerturbation {

    StrainParameterisationPerturbation()=default;

    StrainParameterisationPerturbation(const unsigned int t_number_of_Chebyshev_points) : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
    {}

    StrainParameterisationPerturbation(unsigned int t_ne, const std::vector<bool> &t_admitted_deformations,
                                       const unsigned int t_number_of_Chebyshev_points) :   m_ne(t_ne), m_admitted_deformations(t_admitted_deformations),
                                                                                            m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
    {}


    StrainParameterisationPerturbation(const unsigned int t_number_of_Chebyshev_points,
                                       const std::function<double(const unsigned int, const double&)> t_polynomial_base) :  m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
                                                                                                                           m_polynomial_base(t_polynomial_base)
    {}

    StrainParameterisationPerturbation(unsigned int t_ne, unsigned int t_na,
                                       unsigned int t_number_of_Chebyshev_points,
                                       const std::function<double(const unsigned int, const double&)> t_polynomial_base) :  m_na(t_na), m_ne(t_ne),
                                                                                                                            m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
                                                                                                                            m_polynomial_base(t_polynomial_base)
    {}


    void update(const Eigen::VectorXd &t_Delta_qe,
                const Eigen::VectorXd &t_Delta_dot_qe,
                const Eigen::VectorXd &t_Delta_ddot_qe)
    {

        Eigen::VectorXd Delta_xi;
        Eigen::VectorXd Delta_dot_xi;
        Eigen::VectorXd Delta_ddot_xi;

        Eigen::MatrixXd Delta_xi_stack(6,m_stacks_dimension);

        for(unsigned int i=0; i<m_stacks_dimension; i++){

            Delta_xi = m_B*m_Phi_stack[i]*t_Delta_qe;
            Delta_dot_xi = m_B*m_Phi_stack[i]*t_Delta_dot_qe;
            Delta_ddot_xi = m_B*m_Phi_stack[i]*t_Delta_ddot_qe;

            Delta_xi_stack.col(i) = Delta_xi;

            m_Delta_K_stack->at(i) = Delta_xi.block<3,1>(0,0);
            m_Delta_dot_K_stack->at(i) = Delta_dot_xi.block<3,1>(0,0);
            m_Delta_ddot_K_stack->at(i) = Delta_ddot_xi.block<3,1>(0,0);

            m_Delta_Lambda_stack->at(i) = Delta_xi.block<3,1>(3,0);
            m_Delta_dot_Lambda_stack->at(i) = Delta_dot_xi.block<3,1>(3,0);
            m_Delta_ddot_Lambda_stack->at(i) = Delta_ddot_xi.block<3,1>(3,0);
        }

        //std::cout << "Delta xi : \n" << Delta_xi_stack << "\n\n";
    }





    const std::vector<bool> m_admitted_deformations { false, true, true, false, false, true };

    const unsigned int m_na { [&]()->unsigned int{ return std::count(m_admitted_deformations.begin(), m_admitted_deformations.end(), true);}() };

    const unsigned int m_ne { 3 };

    const Eigen::MatrixXd m_B { [&](){
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);

            std::vector<int> indexes;
            std::for_each(m_admitted_deformations.begin(),
                          m_admitted_deformations.end(),
                          [&indexes, index=0](const bool dof)mutable{   if(dof == true)
                                                                            indexes.push_back(index);
                                                                        index++;});
            const Eigen::MatrixXd map = I(Eigen::all, indexes);
            return map;
                              }() };

    const Eigen::MatrixXd m_B_bar { [&](){

            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);

            std::vector<int> indexes;
            std::for_each(m_admitted_deformations.begin(),
                          m_admitted_deformations.end(),
                          [&indexes, index=0](const bool dof)mutable{   if(dof == false)
                                                                            indexes.push_back(index);
                                                                        index++;});
            const Eigen::MatrixXd map = I(Eigen::all, indexes);
            return map;
                              }() };

    const unsigned int m_number_of_Chebyshev_points { 17 };

    const std::function<double(const unsigned int, const double&)> m_polynomial_base { [](const unsigned int t_point, const double& t_x) {return boost::math::legendre_p(t_point, t_x);} };

    const std::vector<Eigen::MatrixXd> m_Phi_stack { generatePhiStack(m_ne, m_na, ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points), m_polynomial_base) };

    const unsigned int m_stacks_dimension { static_cast<unsigned int>( m_Phi_stack.size() ) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Delta_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Delta_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Delta_dot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Delta_dot_Lambda_stack{ std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };

    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Delta_ddot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Delta_ddot_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_stacks_dimension) };
};



}   //  namespace StrainParameterisation

}   //  namespace CROSP




#endif // STRAIN_PARAMETERISATION_HPP
