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




/*!
 * \brief The PolynomialRepresentation struct contains the reduction of the strain into a polynomial base
 *
 * This struct contains the choosen reduction in order to represent the strain vector field.
 */
struct PolynomialRepresentation {
    PolynomialRepresentation()=default;


    PolynomialRepresentation(const unsigned int t_ne);


    PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations);


    PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations,
                             const unsigned int t_ne);

    PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations,
                             const unsigned int t_ne,
                             const PolynomialBase t_polynomial_base);

    unsigned int getCoordinatesDimension()const{return m_ne*m_na;}



    /*!
     * \brief getPhi this function computes the base matrix Phi
     * \param t_ne the number of modes per admitted deformation
     * \param t_na the number of admitted deformations
     * \param t_X the current abscissa along the normalized arc-length
     * \param t_polynomial_base the polynomial base used to define the matrix
     * \param t_begin the begin of the rod domain
     * \param t_end the end of the rod domain
     * \return the matrix of bases function, namely Phi
     */
    Eigen::MatrixXd getPhi(const double& t_X,
                           const double& t_begin=0,
                           const double& t_end=1)const;







    const std::array<bool, 6> m_admitted_deformations { true, true, true, false, false, false };

    const unsigned int m_na { [&]()->unsigned int{ return std::count(m_admitted_deformations.begin(), m_admitted_deformations.end(), true);}() };

    const unsigned int m_ne { 3 };

    const PolynomialBase m_polynomial_base { legendre_polynomial_base };



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
    }()};


    const Eigen::MatrixXd m_Bbar { [&](){

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






/*!
 * \brief generatePhiStack generates a stack of Phi matrices, one for every Chebyshev points and relative to its abscissa
 * \param t_ne the number of modes per admitted deformation
 * \param t_na the number of admitted deformations
 * \param t_Chebyshev_points the number of Chebyshev points used to discretise the rod
 * \param t_polynomial_base the polynomial base used to define the matrix
 * \return the stack of Phi matrices
 */
std::vector<Eigen::MatrixXd> generatePhiStack(const PolynomialRepresentation &t_polynomial_representation,
                                              const std::vector<double> &t_Chebyshev_points);




}   //  namespace CROSP::polynomial_representation

#endif // POLYNOMIAL_RPRESENTATION_HPPS
