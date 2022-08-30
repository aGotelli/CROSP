#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace CROSP::strain_parameterisation {


/*!
 * \brief StrainParameterisation initialises the object
 * \param t_number_of_Chebyshev_points the number of chebyshev points in which the rod is observed
 */
StrainParameterisation::StrainParameterisation(const unsigned int t_number_of_Chebyshev_points)
    : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


/*!
 * \brief StrainParameterisation initialises the object
 * \param t_ne the number of modes per deformation
 * \param t_number_of_Chebyshev_points the number of chebyshev points in which the rod is observed
 */
StrainParameterisation::StrainParameterisation(const unsigned int t_ne,
                                               const unsigned int t_number_of_Chebyshev_points)
    : m_ne(t_ne),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


/*!
 * \brief StrainParameterisation
 * \param t_constant_strain
 */
StrainParameterisation::StrainParameterisation(const Eigen::VectorXd &t_constant_strain)
    :   m_constrained_strain( defineConstrainedStrain(t_constant_strain) )
{}


/*!
 * \brief StrainParameterisation
 * \param t_ne the number of modes per degree of deformations
 * \param t_admitted_deformations the admitted deformations in the rod
 */
StrainParameterisation::StrainParameterisation(const unsigned int t_ne,
                                               std::array<bool, 6> &t_admitted_deformations)
    : m_ne(t_ne),
      m_admitted_deformations(t_admitted_deformations)
{}

/*!
 * \brief StrainParameterisation
 * \param t_ne
 * \param t_admitted_deformations
 * \param t_number_of_Chebyshev_points
 */
StrainParameterisation::StrainParameterisation(const unsigned int t_ne,
                                               const std::array<bool, 6> &t_admitted_deformations,
                                               const unsigned int t_number_of_Chebyshev_points)
    : m_ne(t_ne),
      m_admitted_deformations(t_admitted_deformations),
      m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}



StrainParameterisation::StrainParameterisation(const unsigned int t_ne,
                                               const std::array<bool, 6> &t_admitted_deformations,
                                               const Eigen::VectorXd &t_constant_strain)
    : m_ne(t_ne),
      m_admitted_deformations(t_admitted_deformations),
      m_constrained_strain( defineConstrainedStrain(t_constant_strain) )
{}


StrainParameterisation::StrainParameterisation(const unsigned int t_ne,
                                               const std::array<bool, 6> &t_admitted_deformations,
                                               const Eigen::VectorXd &t_constant_strain,
                                               const unsigned int t_number_of_Chebyshev_points)
    :   m_ne(t_ne),
        m_admitted_deformations(t_admitted_deformations),
        m_constrained_strain( defineConstrainedStrain(t_constant_strain) ),
        m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
{}


//StrainParameterisation::StrainParameterisation(const unsigned int t_ne,
//                                               const std::array<bool, 6> &t_admitted_deformations,
//                                               const Eigen::VectorXd &t_constant_strain,
//                                               const base_maps::BaseFunction t_polynomial_base)
//    :   m_ne(t_ne),
//        m_admitted_deformations(t_admitted_deformations),
//        m_constrained_strain( defineConstrainedStrain(t_constant_strain) ),
//        m_Phi_stack( base_maps::generatePhiStack(m_ne,
//                                                 m_na,
//                                                 ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points),
//                                                                                     t_polynomial_base) )
//{}


//StrainParameterisation::StrainParameterisation(const unsigned int t_ne,
//                                               const std::array<bool, 6> &t_admitted_deformations,
//                                               const Eigen::VectorXd &t_constant_strain,
//                                               const base_maps::BaseFunction t_polynomial_base,
//                                               const unsigned int t_number_of_Chebyshev_points)
//    :   m_ne(t_ne),
//        m_admitted_deformations(t_admitted_deformations),
//        m_constrained_strain( defineConstrainedStrain(t_constant_strain) ),
//        m_Phi_stack( base_maps::generatePhiStack(m_ne,
//                                                 m_na,
//                                                 ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points),
//                                                                                     t_polynomial_base) ),
//        m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
//{}




void StrainParameterisation::updateStacks(const Eigen::VectorXd &t_qe,
                                          const Eigen::VectorXd &t_dot_qe,
                                          const Eigen::VectorXd &t_ddot_qe)
{

    Eigen::VectorXd xi;
    Eigen::VectorXd dot_xi;
    Eigen::VectorXd ddot_xi;

    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){

        xi = m_map_to_strain_stack[i]*t_qe + m_B_bar*m_constrained_strain;
        dot_xi = m_map_to_strain_stack[i]*t_dot_qe;
        ddot_xi = m_map_to_strain_stack[i]*t_ddot_qe;

        m_K_stack->at(i) = xi.block<3,1>(0,0);
        m_dot_K_stack->at(i) = dot_xi.block<3,1>(0,0);
        m_ddot_K_stack->at(i) = ddot_xi.block<3,1>(0,0);

        m_Lambda_stack->at(i) = xi.block<3,1>(3,0);
        m_dot_Lambda_stack->at(i) = dot_xi.block<3,1>(3,0);
        m_ddot_Lambda_stack->at(i) = ddot_xi.block<3,1>(3,0);

    }

}



Eigen::VectorXd StrainParameterisation::defineConstrainedStrain(const Eigen::VectorXd &t_constant_strain)const
{
    std::vector<int> indexes;
    std::for_each(m_admitted_deformations.begin(),
                  m_admitted_deformations.end(),
                  [&indexes, index=0](const bool dof)mutable{   if(dof == false)
                                                                    indexes.push_back(index);
                                                                index++;});
    Eigen::VectorXd xi_c = t_constant_strain(indexes);

    return xi_c;
}








}   //  namespace CROSP::StrainParameterisation
