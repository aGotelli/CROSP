#ifndef BASE_HPP
#define BASE_HPP

#include <functional>
#include <boost/math/special_functions.hpp>
#include <Eigen/Dense>

namespace CROSP::base {



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



Eigen::MatrixXd getPhi(const unsigned int t_ne,
                              const unsigned int t_na,
                              const double& t_X,
                              const BaseFunction t_polynomial_base=legendre_polynomial_base,
                              const double& t_begin=0,
                              const double& t_end=1);

std::vector<Eigen::MatrixXd> generatePhiStack(const unsigned int t_ne,
                                              const unsigned int t_na,
                                              const std::vector<double> &t_Chebyshev_points,
                                              const BaseFunction t_polynomial_base=legendre_polynomial_base);


}   //  namespace CROSP::base

#endif // BASE_HPP
