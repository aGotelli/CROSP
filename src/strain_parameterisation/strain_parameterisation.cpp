#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"

#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace CROSP::strain_parameterisation {


StrainParameterisation::StrainParameterisation(const Eigen::VectorXd t_constant_strain)
    : m_constant_strain(t_constant_strain)
{}



}   //  namespace CROSP::StrainParameterisation
