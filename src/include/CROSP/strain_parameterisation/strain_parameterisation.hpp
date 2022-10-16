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

#include "CROSP/polynomial_representation/polynomial_representation.hpp"

namespace CROSP::strain_parameterisation {






/*!
 * \brief The StrainParameterisation class contains the strains and its derivatives with respect to time
 *
 * This class contains all the terms related to the rod strain, with the function needed to update such strain.
 * It contains the strains with their first and second derivatives.
 *
 * This class is build on the assumption that the rod is observed with a fixed number of Chebyshev points.
 * As a result, all the strain and their derivatives are evaluated at these Chebyshev points and stored in
 * a std::vector.
 */
class StrainParameterisation {

public:

    StrainParameterisation()=default;


    StrainParameterisation(const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(const Eigen::VectorXd &t_constant_strain);


    StrainParameterisation(const polynomial_representation::PolynomialRepresentation t_polynomial_representation);


    StrainParameterisation(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                           const unsigned int t_number_of_Chebyshev_points);

    StrainParameterisation(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                           const Eigen::VectorXd &t_constant_strain);


    StrainParameterisation(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                           const Eigen::VectorXd &t_constant_strain,
                           const unsigned int t_number_of_Chebyshev_points);




    ~StrainParameterisation()=default;


    /*!
     * \brief updateStacks updates the stacks of strain and its first and second derivative
     * \param t_qe the generalised elastic coordinates
     * \param t_dot_qe the first derivative of the generalised elastic coordinates
     * \param t_ddot_qe the second derivative of the generalised elastic coordinates
     */
    void updateStacks(const Eigen::VectorXd &t_qe,
                      const Eigen::VectorXd &t_dot_qe,
                      const Eigen::VectorXd &t_ddot_qe);



    //  The number of Chebyshev points used to represent the rod
    const unsigned int m_number_of_Chebyshev_points { 17 };


    //  The stacks for the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };

    //  The stacks for the first derivative of the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Lambda_stack{ std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };

    //  The stacks for the first derivative of the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack { std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points) };

private:

    //  The polynomial representation of the field of strain
    const polynomial_representation::PolynomialRepresentation m_polynomial_representation { polynomial_representation::PolynomialRepresentation() };


public:

    //  The stack of Phi along the Chebyshev points
    const std::vector<Eigen::MatrixXd> m_Phi_stack { polynomial_representation::generatePhiStack(m_polynomial_representation,
                                                                                                 ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points)) };


    //  The vector stack of B Phi used to map the generalised coordinates into the strain field
    const std::vector<Eigen::MatrixXd> m_map_to_strain_stack { [&](){
            std::vector<Eigen::MatrixXd> map_to_strain_stack(m_number_of_Chebyshev_points);

            const std::vector<Eigen::MatrixXd> Phi_stack =
                    polynomial_representation::generatePhiStack(m_polynomial_representation,
                                                                ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points));


            std::generate(map_to_strain_stack.begin(), map_to_strain_stack.end(), [&, index=0]()mutable{
               return m_polynomial_representation.m_B*Phi_stack[index++];
            });

            return map_to_strain_stack;
        }() };
private:

    /*!
     * \brief defineConstrainedStrain defines the constrained strains without the degrees of freedom
     * \param t_constant_strain the constant strain along the rod with is 6 components
     * \return the constrained strain of dimension 6-na
     *
     * This function takes the constant strain as a full 6x1 vector and returns the corresponding
     * subvector of dimension (6-na)x1 that contains the components that do not belong to the
     * rod allowed deformations.
     */
    Eigen::VectorXd defineConstrainedStrain(const Eigen::VectorXd &t_constant_strain)const;


    //  The constrained strain xi_c
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


};




}   //  namespace CROSP::StrainParameterisation





#endif // STRAIN_PARAMETERISATION_HPP
