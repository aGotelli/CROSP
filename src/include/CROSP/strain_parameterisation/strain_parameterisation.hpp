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

#include <boost/math/special_functions/legendre.hpp>

#include "math_tools/Chebyshev/chebyshev_differentiation.hpp"
#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"



namespace CROSP::strain_parameterisation {


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





struct StrainParameterisationBase {

    StrainParameterisationBase()=default;


    StrainParameterisationBase(const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisationBase(unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisationBase(unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const unsigned int t_number_of_Chebyshev_points,
                           const Eigen::VectorXd &t_constrained_strain,
                           const BaseFunction t_polynomial_base=default_base_function);

    ~StrainParameterisationBase()=default;


    inline unsigned int getNumberOfChebyshewPoints()const{return m_number_of_Chebyshev_points;}

    inline unsigned int getCoordinatesDimention()const{return m_na*m_ne;}


    void updateStacks(const Eigen::VectorXd &t_qe,
                      const Eigen::VectorXd &t_dot_qe,
                      const Eigen::VectorXd &t_ddot_qe);

    virtual void updateStrain(const Eigen::VectorXd &t_qe,
                              const Eigen::VectorXd &t_dot_qe,
                              const ::LieAlgebra::Vector6d &t_xi,
                              const Eigen::VectorXd &t_ddot_qe,
                              const ::LieAlgebra::Vector6d &t_dot_xi,
                              const ::LieAlgebra::Vector6d &t_ddot_xi)const=0;








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

struct StrainParameterisation {

    StrainParameterisation()=default;


    StrainParameterisation(const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(unsigned int t_ne,
                           const std::vector<bool> &t_admitted_deformations,
                           const unsigned int t_number_of_Chebyshev_points,
                           const Eigen::VectorXd &t_constrained_strain,
                           const BaseFunction t_polynomial_base=default_base_function);

    ~StrainParameterisation()=default;


    inline unsigned int getNumberOfChebyshewPoints()const{return m_number_of_Chebyshev_points;}

    inline unsigned int getCoordinatesDimention()const{return m_na*m_ne;}


    virtual void update(const Eigen::VectorXd &t_qe,
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


struct StrainParameterisationDelta : public StrainParameterisation {

    StrainParameterisationDelta()=default;

    StrainParameterisationDelta(const unsigned int t_number_of_Chebyshev_points);

    StrainParameterisationDelta(unsigned int t_ne,
                                const std::vector<bool> &t_admitted_deformations,
                                const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisationDelta(unsigned int t_ne,
                                const std::vector<bool> &t_admitted_deformations,
                                const unsigned int t_number_of_Chebyshev_points,
                                const Eigen::VectorXd &t_constrained_strain,
                                const BaseFunction t_polynomial_base=default_base_function);

    virtual void update(const Eigen::VectorXd &t_Delta_qe,
                        const Eigen::VectorXd &t_Delta_dot_qe,
                        const Eigen::VectorXd &t_Delta_ddot_qe)override;

};




}   //  namespace CROSP::StrainParameterisation





#endif // STRAIN_PARAMETERISATION_HPP
