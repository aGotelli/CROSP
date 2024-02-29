/**
 * \file polynomial_representation.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the definition of class and function used for the parameterisation of the strain field with a polynomial basis
 * \date 12-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */

#ifndef POLYNOMIAL_RPRESENTATION_HPP
#define POLYNOMIAL_RPRESENTATION_HPP

#include <functional>
#include <boost/math/special_functions.hpp>
#include <Eigen/Dense>

#include <numeric>
#include <memory>

namespace CROSP::polynomial_representation {


typedef std::function<double(const unsigned int, const double&)> PolynomialBase;


/// legendre_polynomial_base defines a base function using the Legendre polynomials
static const PolynomialBase legendre_polynomial_base { [](const unsigned int t_point,const double& t_x)
                                                        {
                                                            return boost::math::legendre_p(t_point, t_x);
                                                        }
                                                    };

/// chebyshev_polynomial_base defines a base funtion using the chebyshev polynomials
static const PolynomialBase chebyshev_polynomial_base { [](const unsigned int t_point,const double& t_x)
                                                        {
                                                            return boost::math::chebyshev_t(t_point, t_x);
                                                        }
                                                    };


static constexpr unsigned int default_number_of_modes = 4;




/*!
 * \brief The PolynomialRepresentation struct contains the reduction of the strain into a polynomial base
 *
 * This struct contains the choosen reduction in order to represent the strain vector field.
 */
struct PolynomialRepresentation {

    PolynomialRepresentation()=default;


    PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations,
                             const unsigned int t_number_of_modes,
                             const PolynomialBase t_polynomial_base=legendre_polynomial_base);


    PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations,
                             const std::vector<unsigned int> t_number_of_modes_stack,
                             const PolynomialBase t_polynomial_base=legendre_polynomial_base);

    PolynomialRepresentation(const PolynomialRepresentation &t_other);


    /*!
     * \brief getCoordinatesDimension computes the dimention of q_e the vector of generalised coordinates describing thr rod shape
     * \return the dimention of q_e the vector of generalised coordinates describing thr rod shape
     */
    unsigned int getCoordinatesDimension()const{return m_total_number_of_modes;}



    /*!
     * \brief getPhi this function computes the base matrix Phi
     * \param t_X the current abscissa along the normalized arc-length
     * \param t_begin the begin of the rod domain (default 0)
     * \param t_end the end of the rod domain (default 1)
     * \return the matrix of bases function, namely Phi
     */
    Eigen::MatrixXd getPhi(const double& t_X,
                           const double& t_begin=0,
                           const double& t_end=1)const;






    //  A vector representing the allowed deformation of the rod (default kirkhoff rod)
    std::array<bool, 6> m_admitted_deformations { false, false, true, false, false, false };

    //  The number of allowed deformation (asserted from the vector)
    unsigned int m_na { [&]()->unsigned int{ return std::count(m_admitted_deformations.begin(),
                                                                     m_admitted_deformations.end(), true);}() };

    //  Number of elastic modes per allowed deformation
    std::vector<unsigned int> m_number_of_modes_stack { default_number_of_modes };
    unsigned int m_total_number_of_modes { static_cast<unsigned int>(std::accumulate(m_number_of_modes_stack.begin(),
                                                                                           m_number_of_modes_stack.end(), 0)) };

    //  The polynomial base used to discretize the strain field
    PolynomialBase m_polynomial_base { legendre_polynomial_base };

    //  The matrix mapping the allowed strain into the full strain field
    Eigen::MatrixXd m_B { [&](){

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);

        std::vector<int> indexes;
        std::for_each(m_admitted_deformations.begin(),
                      m_admitted_deformations.end(),
                      [&indexes, index=0](const bool dof)mutable{   if(dof == true)
                                                                        indexes.push_back(index);
                                                                    index++;});
        const Eigen::MatrixXd map = I(Eigen::all, indexes);
        return map;
    }()};


    //  The matrix mapping the constrained strain into the full strain field
    Eigen::MatrixXd m_Bbar { [&](){

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);

        std::vector<int> indexes;
        std::for_each(m_admitted_deformations.begin(),
                      m_admitted_deformations.end(),
                      [&indexes, index=0](const bool dof)mutable{   if(dof == false)
                                                                        indexes.push_back(index);
                                                                    index++;});
        const Eigen::MatrixXd map = I(Eigen::all, indexes);
        return map;
    }()};




};






}   //  namespace CROSP::polynomial_representation

#endif // POLYNOMIAL_RPRESENTATION_HPP
