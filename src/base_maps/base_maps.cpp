#include "CROSP/base_maps/base_maps.hpp"

#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace CROSP::base_maps {



Eigen::MatrixXd getPhi(const unsigned int t_ne,
                       const unsigned int t_na,
                       const double& t_X,
                       const PolynomialBase t_polynomial_base,
                       const double& t_begin,
                       const double& t_end)
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

std::vector<Eigen::MatrixXd> generatePhiStack(const unsigned int t_ne,
                                              const unsigned int t_na,
                                              const std::vector<double> &t_Chebyshev_points,
                                              const PolynomialBase t_polynomial_base)
{
    std::vector<Eigen::MatrixXd> Phi_stack( t_Chebyshev_points.size() );

    std::generate(Phi_stack.begin(), Phi_stack.end(), [&, index=0]()mutable{
        const auto Phi = getPhi(t_ne, t_na, t_Chebyshev_points[index], t_polynomial_base);
        index++;
        return Phi;
    });

    return Phi_stack;
}


Eigen::MatrixXd getB(const std::array<bool, 6> &t_allowed_deformations)
{
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);

    std::vector<int> indexes;
    std::for_each(t_allowed_deformations.begin(),
                  t_allowed_deformations.end(),
                  [&indexes, index=0](const bool dof)mutable{   if(dof == true)
                                                                    indexes.push_back(index);
                                                                index++;});
    const Eigen::MatrixXd map = I(Eigen::all, indexes);
    return map;
}


Eigen::MatrixXd getBbar(const std::array<bool, 6> &t_allowed_deformations)
{
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);

    std::vector<int> indexes;
    std::for_each(t_allowed_deformations.begin(),
                  t_allowed_deformations.end(),
                  [&indexes, index=0](const bool dof)mutable{   if(dof == false)
                                                                    indexes.push_back(index);
                                                                index++;});
    const Eigen::MatrixXd map = I(Eigen::all, indexes);
    return map;
}


}   //  namespace CROSP::base_maps
