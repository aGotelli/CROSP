#ifndef PARAMETERISATIONSTACK_HPP
#define PARAMETERISATIONSTACK_HPP

#include <vector>
#include <Eigen/Dense>

#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/polynomial_representation/polynomial_representation.hpp"


namespace CROSP::spectral_integrators {


struct ParameterisationStack {


    ParameterisationStack(const strain_parameterisation::StrainParameterisation &t_strain_parameterisation,
                          const polynomial_representation::PolynomialRepresentation &t_polynomial_representation,
                          const unsigned int t_number_of_Chebyshev_points);


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


    //  The strain parameterisation of the rod
    strain_parameterisation::StrainParameterisation m_strain_parameterisation;

    //  The polynomial representation of the field of strain
    polynomial_representation::PolynomialRepresentation m_polynomial_representation;

    //  The stacks for the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Gamma_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };

    //  The stacks for the first derivative of the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Gamma_stack{
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };

    //  The stacks for the first derivative of the strain, decomposed in angular and linear part
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Gamma_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>(m_number_of_Chebyshev_points)
    };







    //  The vector stack of B Phi used to map the generalised coordinates into the strain field
    std::vector<Eigen::MatrixXd> m_map_to_strain_stack { [this](){

        std::vector<Eigen::MatrixXd> map_to_strain_stack(m_number_of_Chebyshev_points);

        std::vector<Eigen::MatrixXd> Phi_stack( m_number_of_Chebyshev_points );

        const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints( m_number_of_Chebyshev_points );

        std::generate(Phi_stack.begin(), Phi_stack.end(), [&, index=0]()mutable{
            const auto Phi = m_polynomial_representation.getPhi( Chebyshev_points[index] );
            index++;
            return Phi;
        });


        std::generate(map_to_strain_stack.begin(), map_to_strain_stack.end(), [&, index=0]()mutable{
            return m_polynomial_representation.m_B*Phi_stack[index++];
        });

        return map_to_strain_stack;
    }() };


};

}   //  namespace   CROSP::spectral_integrators



#endif // PARAMETERISATIONSTACK_HPP
