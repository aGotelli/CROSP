#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace CROSP::strain_parameterisation {


Eigen::MatrixXd getPhi(const unsigned int t_ne,
                       const unsigned int t_na,
                       const double& t_X,
                       const BaseFunction t_polynomial_base,
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
                                              const BaseFunction t_polynomial_base)
{
    std::vector<Eigen::MatrixXd> Phi_stack( t_Chebyshev_points.size() );

    std::generate(Phi_stack.begin(), Phi_stack.end(), [&, index=0]()mutable{
        const auto Phi = getPhi(t_ne, t_na, t_Chebyshev_points[index], t_polynomial_base);
        index++;
        return Phi;
    });

    return Phi_stack;
}



StrainParameterisationBase::StrainParameterisationBase(const unsigned int t_number_of_Chebyshev_points) :
    m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


StrainParameterisationBase::StrainParameterisationBase(unsigned int t_ne,
                                               const std::vector<bool> &t_admitted_deformations,
                                               const unsigned int t_number_of_Chebyshev_points) :
    m_ne(t_ne),
    m_admitted_deformations(t_admitted_deformations),
    m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


StrainParameterisationBase::StrainParameterisationBase(unsigned int t_ne,
                                               const std::vector<bool> &t_admitted_deformations,
                                               const unsigned int t_number_of_Chebyshev_points,
                                               const Eigen::VectorXd &t_constrained_strain,
                                               const std::function<double(const unsigned int, const double&)> t_polynomial_base) :
    m_ne(t_ne),
    m_admitted_deformations(t_admitted_deformations),
    m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
    m_constrained_strain(t_constrained_strain),
    m_polynomial_base(t_polynomial_base)
{}




void StrainParameterisationBase::updateStacks(const Eigen::VectorXd &t_qe,
            const Eigen::VectorXd &t_dot_qe,
            const Eigen::VectorXd &t_ddot_qe)
{

    Eigen::VectorXd xi;
    Eigen::VectorXd dot_xi;
    Eigen::VectorXd ddot_xi;

    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){

        xi = m_strains_map_stack[i]*t_qe + m_B_bar*m_constrained_strain;
        dot_xi = m_strains_map_stack[i]*t_dot_qe;
        ddot_xi = m_strains_map_stack[i]*t_ddot_qe;

        m_K_stack->at(i) = xi.block<3,1>(0,0);
        m_dot_K_stack->at(i) = dot_xi.block<3,1>(0,0);
        m_ddot_K_stack->at(i) = ddot_xi.block<3,1>(0,0);

        m_Lambda_stack->at(i) = xi.block<3,1>(3,0);
        m_dot_Lambda_stack->at(i) = dot_xi.block<3,1>(3,0);
        m_ddot_Lambda_stack->at(i) = ddot_xi.block<3,1>(3,0);

    }

}



StrainParameterisation::StrainParameterisation(const unsigned int t_number_of_Chebyshev_points) :
    m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


StrainParameterisation::StrainParameterisation(unsigned int t_ne,
                                               const std::vector<bool> &t_admitted_deformations,
                                               const unsigned int t_number_of_Chebyshev_points) :
    m_ne(t_ne),
    m_admitted_deformations(t_admitted_deformations),
    m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


StrainParameterisation::StrainParameterisation(unsigned int t_ne,
                                               const std::vector<bool> &t_admitted_deformations,
                                               const unsigned int t_number_of_Chebyshev_points,
                                               const Eigen::VectorXd &t_constrained_strain,
                                               const std::function<double(const unsigned int, const double&)> t_polynomial_base) :
    m_ne(t_ne),
    m_admitted_deformations(t_admitted_deformations),
    m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
    m_constrained_strain(t_constrained_strain),
    m_polynomial_base(t_polynomial_base)
{}




void StrainParameterisation::update(const Eigen::VectorXd &t_qe,
            const Eigen::VectorXd &t_dot_qe,
            const Eigen::VectorXd &t_ddot_qe)
{

    Eigen::VectorXd xi;
    Eigen::VectorXd dot_xi;
    Eigen::VectorXd ddot_xi;

    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){

        xi = m_strains_map_stack[i]*t_qe + m_B_bar*m_constrained_strain;
        dot_xi = m_strains_map_stack[i]*t_dot_qe;
        ddot_xi = m_strains_map_stack[i]*t_ddot_qe;

        m_K_stack->at(i) = xi.block<3,1>(0,0);
        m_dot_K_stack->at(i) = dot_xi.block<3,1>(0,0);
        m_ddot_K_stack->at(i) = ddot_xi.block<3,1>(0,0);

        m_Lambda_stack->at(i) = xi.block<3,1>(3,0);
        m_dot_Lambda_stack->at(i) = dot_xi.block<3,1>(3,0);
        m_ddot_Lambda_stack->at(i) = ddot_xi.block<3,1>(3,0);

    }

}



StrainParameterisationDelta::StrainParameterisationDelta(const unsigned int t_number_of_Chebyshev_points) :
    StrainParameterisation(t_number_of_Chebyshev_points)
{}

StrainParameterisationDelta::StrainParameterisationDelta(unsigned int t_ne,
                                                         const std::vector<bool> &t_admitted_deformations,
                                                         const unsigned int t_number_of_Chebyshev_points) :
    StrainParameterisation(t_ne,
                           t_admitted_deformations,
                           t_number_of_Chebyshev_points)
{}



StrainParameterisationDelta::StrainParameterisationDelta(unsigned int t_ne,
                                const std::vector<bool> &t_admitted_deformations,
                                const unsigned int t_number_of_Chebyshev_points,
                                const Eigen::VectorXd &t_constrained_strain,
                                const BaseFunction t_polynomial_base) :
    StrainParameterisation(t_ne,
                           t_admitted_deformations,
                           t_number_of_Chebyshev_points,
                           t_constrained_strain,
                           t_polynomial_base)
{}

void StrainParameterisationDelta::update(const Eigen::VectorXd &t_Delta_qe,
                                         const Eigen::VectorXd &t_Delta_dot_qe,
                                         const Eigen::VectorXd &t_Delta_ddot_qe)
{
    Eigen::VectorXd Delta_xi;
    Eigen::VectorXd Delta_dot_xi;
    Eigen::VectorXd Delta_ddot_xi;


    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){

        Delta_xi = m_B*m_Phi_stack[i]*t_Delta_qe;
        Delta_dot_xi = m_B*m_Phi_stack[i]*t_Delta_dot_qe;
        Delta_ddot_xi = m_B*m_Phi_stack[i]*t_Delta_ddot_qe;


        m_K_stack->at(i) = Delta_xi.block<3,1>(0,0);
        m_dot_K_stack->at(i) = Delta_dot_xi.block<3,1>(0,0);
        m_ddot_K_stack->at(i) = Delta_ddot_xi.block<3,1>(0,0);

        m_Lambda_stack->at(i) = Delta_xi.block<3,1>(3,0);
        m_dot_Lambda_stack->at(i) = Delta_dot_xi.block<3,1>(3,0);
        m_ddot_Lambda_stack->at(i) = Delta_ddot_xi.block<3,1>(3,0);
    }
}





}   //  namespace CROSP::StrainParameterisation
