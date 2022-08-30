#include "CROSP/base_maps/base_maps.hpp"

#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace CROSP::base_maps {





PolynomialRepresentation::PolynomialRepresentation(const unsigned int t_ne)
    : m_ne(t_ne)
{}



PolynomialRepresentation::PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations)
    : m_admitted_deformations(t_admitted_deformations)
{}



PolynomialRepresentation::PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations,
                                                   const unsigned int t_ne)
    : m_admitted_deformations(t_admitted_deformations),
      m_ne(t_ne)
{}


PolynomialRepresentation::PolynomialRepresentation(const std::array<bool, 6> t_admitted_deformations,
                                                   const unsigned int t_ne,
                                                   const PolynomialBase t_polynomial_base)
    : m_admitted_deformations(t_admitted_deformations),
      m_ne(t_ne),
      m_polynomial_base(t_polynomial_base)
{}




Eigen::MatrixXd PolynomialRepresentation::getPhi(const double& t_X,
                                                 const double& t_begin,
                                                 const double& t_end) const
{
    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    double x = ( 2 * t_X - ( t_end + t_begin) ) / ( t_end - t_begin );

    //  Compute the values of the polynomial for every element of the strain field
    Eigen::VectorXd Phi_i(m_ne, 1);
    for(unsigned int i=0; i<m_ne; i++)
        Phi_i[i] = m_polynomial_base(i, x);


    //  Define the matrix of bases
    Eigen::MatrixXd Phi = Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_na, m_na), Phi_i.transpose());


    return Phi;
}






std::vector<Eigen::MatrixXd> generatePhiStack(const PolynomialRepresentation &t_polynomial_representation,
                                              const std::vector<double> &t_Chebyshev_points)
{
    std::vector<Eigen::MatrixXd> Phi_stack( t_Chebyshev_points.size() );

    std::generate(Phi_stack.begin(), Phi_stack.end(), [&, index=0]()mutable{
        const auto Phi = t_polynomial_representation.getPhi(t_Chebyshev_points[index]);
        index++;
        return Phi;
    });

    return Phi_stack;
}



}   //  namespace CROSP::base_maps
