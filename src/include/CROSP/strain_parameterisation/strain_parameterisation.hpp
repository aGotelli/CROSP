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


#include "math_tools/Chebyshev/chebyshev_differentiation.hpp"
#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/base/base.hpp"


namespace CROSP::strain_parameterisation {







class StrainParameterisation {

public:

    StrainParameterisation()=default;


    StrainParameterisation(const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(const unsigned int t_ne,
                           const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(const unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations);


    StrainParameterisation(const unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(const unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const Eigen::VectorXd &t_constant_strain);


    StrainParameterisation(const unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const Eigen::VectorXd &t_constant_strain,
                           const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(const unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const Eigen::VectorXd &t_constant_strain,
                           const base::BaseFunction t_polynomial_base);


    StrainParameterisation(const unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const Eigen::VectorXd &t_constant_strain,
                           const base::BaseFunction t_polynomial_base,
                           const unsigned int t_number_of_Chebyshev_points);

    ~StrainParameterisation()=default;


    inline unsigned int getNumberOfChebyshewPoints()const{return m_number_of_Chebyshev_points;}

    inline unsigned int getCoordinatesDimention()const{return m_na*m_ne;}


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
            Eigen::VectorXd constant_strain(6);
            constant_strain <<   0,
                                 0,
                                 0,
                                 1.0,
                                 0,
                                 0;

            return defineConstrainedStrain(constant_strain);
        }() };



    const unsigned int m_number_of_Chebyshev_points { 17 };

    const std::function<double(const unsigned int, const double&)> m_polynomial_base { [](const unsigned int t_point, const double& t_x) {return boost::math::legendre_p(t_point, t_x);} };

    const std::vector<Eigen::MatrixXd> m_Phi_stack { base::generatePhiStack(m_ne, m_na, ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points), m_polynomial_base) };

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

private:

    Eigen::VectorXd defineConstrainedStrain(const Eigen::VectorXd &t_constant_strain)const;

};




}   //  namespace CROSP::StrainParameterisation





#endif // STRAIN_PARAMETERISATION_HPP
