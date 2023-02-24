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


    StrainParameterisation(const std::shared_ptr<const polynomial_representation::PolynomialRepresentation> t_polynomial_representation,
                           const unsigned int t_number_of_Chebyshev_points);


    StrainParameterisation(const Eigen::VectorXd t_constant_strain,
                           const std::shared_ptr<const strain_parameterisation::StrainParameterisation> t_strain_parameterisation);



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


    const Eigen::VectorXd m_constant_strain {
        (Eigen::VectorXd(6) << 0, 0, 0, 1, 0, 0 ).finished()
    };



    //  The stacks for the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Lambda_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };

    //  The stacks for the first derivative of the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Lambda_stack{
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };

    //  The stacks for the first derivative of the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Lambda_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };



    //  The polynomial representation of the field of strain
    const std::shared_ptr<const polynomial_representation::PolynomialRepresentation> m_polynomial_representation {
        std::make_shared<const polynomial_representation::PolynomialRepresentation>()
    };




//    //  The stack of Phi along the Chebyshev points
//    std::vector<Eigen::MatrixXd> m_Phi_stack {[this](){
//            std::vector<Eigen::MatrixXd> Phi_stack( m_number_of_Chebyshev_points );

//            const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints( m_number_of_Chebyshev_points );

//            std::generate(Phi_stack.begin(), Phi_stack.end(), [&, index=0]()mutable{
//                const auto Phi = m_polynomial_representation->getPhi( Chebyshev_points[index] );
//                index++;
//                return Phi;
//            });

//            return Phi_stack;
//        }() };


    //  The vector stack of B Phi used to map the generalised coordinates into the strain field
    std::vector<Eigen::MatrixXd> m_map_to_strain_stack { [this](){
        std::vector<Eigen::MatrixXd> map_to_strain_stack(m_number_of_Chebyshev_points);

        std::vector<Eigen::MatrixXd> Phi_stack( m_number_of_Chebyshev_points );

        const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints( m_number_of_Chebyshev_points );

        std::generate(Phi_stack.begin(), Phi_stack.end(), [&, index=0]()mutable{
            const auto Phi = m_polynomial_representation->getPhi( Chebyshev_points[index] );
            index++;
            return Phi;
        });


        std::generate(map_to_strain_stack.begin(), map_to_strain_stack.end(), [&, index=0]()mutable{
            return m_polynomial_representation->m_B*Phi_stack[index++];
        });

        return map_to_strain_stack;
    }() };

private:


    //  The constrained strain xi_c
//    const Eigen::VectorXd constant_strain { [&](){

//            //  constrain only the evolution along x
//            Eigen::VectorXd constant_strain(6);
//            constant_strain <<   0,
//                                 0,
//                                 0,
//                                 1.0,
//                                 0,
//                                 0;

//            //  Get the matrix B bar
//            const auto B_bar = m_polynomial_representation->m_Bbar;

//            //  Define the constrained strain from the rod DoFs
//            std::vector<int> indexes;
//            std::for_each(m_polynomial_representation->m_admitted_deformations.begin(),
//                          m_polynomial_representation->m_admitted_deformations.end(),
//                          [&indexes, index=0](const bool dof)mutable{   if(dof == false)
//                                                                            indexes.push_back(index);
//                                                                        index++;});
//            Eigen::VectorXd xi_c = constant_strain(indexes);

//            //  Always check that Gamma x is 1
//            if(xi_c(3) != 1)
//                xi_c(3) = 1.0;


//            return xi_c;
//        }() };


};




}   //  namespace CROSP::StrainParameterisation





#endif // STRAIN_PARAMETERISATION_HPP
