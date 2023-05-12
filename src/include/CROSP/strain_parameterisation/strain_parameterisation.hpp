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
struct StrainParameterisation {


    StrainParameterisation()=default;


    StrainParameterisation(const Eigen::VectorXd t_constant_strain);



    ~StrainParameterisation()=default;



    const Eigen::VectorXd m_constant_strain {
        (Eigen::VectorXd(6) << 0, 0, 0, 1, 0, 0 ).finished()
    };



};




}   //  namespace CROSP::StrainParameterisation





#endif // STRAIN_PARAMETERISATION_HPP
