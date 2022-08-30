#ifndef BASE_MAPS_HPP
#define BASE_MAPS_HPP

#include <functional>
#include <boost/math/special_functions.hpp>
#include <Eigen/Dense>

namespace CROSP::base_maps {



/*!
 * \brief The StrainDiscretization struct contains the reduction of the strain into a polynomial base
 *
 * This struct contains the choosen reduction in order to represent the strain vector field.
 */
struct StrainDiscretization {
    StrainDiscretization()=default;


    StrainDiscretization(const unsigned int t_ne);


    StrainDiscretization(const std::array<bool, 6> t_admitted_deformations);


    StrainDiscretization(const std::array<bool, 6> t_admitted_deformations,
                         const unsigned int t_ne);



    const std::array<bool, 6> m_admitted_deformations { true, true, true, false, false, false };

    const unsigned int m_na { [&]()->unsigned int{ return std::count(m_admitted_deformations.begin(), m_admitted_deformations.end(), true);}() };

    const unsigned int m_ne { 3 };
};



typedef std::function<double(const unsigned int, const double&)> BaseFunction;

/// legendre_polynomial_base defines a base function using the Legendre polynomials
static const BaseFunction legendre_polynomial_base { [](const unsigned int t_point,const double& t_x)
                                                        {
                                                            return boost::math::legendre_p(t_point, t_x);
                                                        }
                                                    };

/// chebyshev_polynomial_base defines a base funtion using the chebyshev polynomials
static const BaseFunction chebyshev_polynomial_base { [](const unsigned int t_point,const double& t_x)
                                                        {
                                                            return boost::math::chebyshev_t(t_point, t_x);
                                                        }
                                                    };


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
Eigen::MatrixXd getPhi(const unsigned int t_ne,
                              const unsigned int t_na,
                              const double& t_X,
                              const BaseFunction t_polynomial_base=legendre_polynomial_base,
                              const double& t_begin=0,
                              const double& t_end=1);


/*!
 * \brief generatePhiStack generates a stack of Phi matrices, one for every Chebyshev points and relative to its abscissa
 * \param t_ne the number of modes per admitted deformation
 * \param t_na the number of admitted deformations
 * \param t_Chebyshev_points the number of Chebyshev points used to discretise the rod
 * \param t_polynomial_base the polynomial base used to define the matrix
 * \return the stack of Phi matrices
 */
std::vector<Eigen::MatrixXd> generatePhiStack(const unsigned int t_ne,
                                              const unsigned int t_na,
                                              const std::vector<double> &t_Chebyshev_points,
                                              const BaseFunction t_polynomial_base=legendre_polynomial_base);


/*!
 * \brief getB computes the map matrix used to map the allowed strain into the full strain field
 * \param t_allowed_deformations the set of allowed deformation
 * \return the map matrix used to map the allowed strain into the full strain field
 */
Eigen::MatrixXd getB(const std::array<bool, 6> &t_allowed_deformations);


/*!
 * \brief getBbar computes the map matrix used to map the constrained strain into the full strain field
 * \param t_allowed_deformations the set of allowed deformation
 * \return the map matrix used to map the constrained strain into the full strain field
 */
Eigen::MatrixXd getBbar(const std::array<bool, 6> &t_allowed_deformations);


}   //  namespace CROSP::base_maps

#endif // BASE_MAPS_HPP
