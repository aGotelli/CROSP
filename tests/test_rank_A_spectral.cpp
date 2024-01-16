

#include "CROSP/polynomial_representation/polynomial_representation.hpp"
#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"





#include "math_tools/Chebyshev/chebyshev_differentiation.hpp"
#include "utilities/Eigen/eigen_io.hpp"


int main(int argc, char *argv[])
{

    const unsigned int number_of_Chebyshev_points = 7;


    const unsigned int ne = 5;


    const std::array<bool, 6> admitted_deformations = {
        true,
        true,
        true,

        false,
        false,
        false
    };

    const unsigned int na = std::count(admitted_deformations.begin(),
                                       admitted_deformations.end(),
                                       true);

    const unsigned int coordinated_dimension = na * ne;





    ::CROSP::polynomial_representation::PolynomialRepresentation
        polynomial_representation(admitted_deformations,
                                  ne,
                                  //::CROSP::polynomial_representation::chebyshev_polynomial_base
                                  ::CROSP::polynomial_representation::legendre_polynomial_base);

    double radius = 0.01;
    ::CROSP::rod_properties::CircularCrossSectionUPtr circular_cross_section =
        std::make_unique<::CROSP::rod_properties::CircularCrossSection>(radius);

    double width = 0.01;
    double height = 0.02;
    ::CROSP::rod_properties::RectangularCrossSectionUPtr rectangular_cross_section =
        std::make_unique<::CROSP::rod_properties::RectangularCrossSection>(width, height);


    double length = 2.45;
    ::CROSP::rod_properties::RodDimensions rod_dimensions(std::move(rectangular_cross_section),
                                                          length);


    ::CROSP::rod_properties::MaterialProperties rod_material{.m_G=0};
    rod_material.m_E = 741e6;


    ::CROSP::strain_parameterisation_stack::StrainParameterisationStack m_strain_parameterisation_stack(polynomial_representation,
                                                                                              number_of_Chebyshev_points);



    Eigen::VectorXd q = Eigen::VectorXd::Random(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(coordinated_dimension);

    // q(0) = 1;
    // q(ne) = 2;
    // q(2*ne) = 3;




    m_strain_parameterisation_stack.updateStrainParameterisation(q, dot_q, ddot_q);
















    unsigned int state_dimension = 3;


    ::Chebyshev::INTEGRATION_DIRECTION integration_direction = ::Chebyshev::INTEGRATION_DIRECTION::FORWARD;

    double upper_integration_limit = 1.0;


    /// \brief m_D_NN is the section of the Differentiation matrix for the unknown points
    Eigen::MatrixXd m_D_NN { ::Chebyshev::getD_NN(number_of_Chebyshev_points, state_dimension, integration_direction) };

    /// \brief m_D_IN is the section of the Differentiation matrix giving the influence of the initial conditions onto the other points
    Eigen::MatrixXd m_D_IN { ::Chebyshev::getD_IN(number_of_Chebyshev_points, state_dimension, integration_direction) };

    /// \brief m_C_NN is the matrix that has to be inverted to solve the system \f$ y = C^{-1} b \f$
    Eigen::MatrixXd m_C_NN { m_D_NN/upper_integration_limit };

    Eigen::MatrixXd A_NN { 0*m_D_NN };


    std::vector<unsigned int> integration_points { ::Chebyshev::defineIntegrationPoints(number_of_Chebyshev_points, integration_direction) };



    unsigned int index = 0;
    for(const auto current_chebyshev_point : integration_points) {
        //  Compute the A matrix of Q' = 1/2 A(K) Q
        Eigen::MatrixXd A_at_Chebyshev_point = ::LieAlgebra::skew(m_strain_parameterisation_stack.m_K_stack->at(current_chebyshev_point));

        for (unsigned int row = 0; row < state_dimension; ++row) {
            for (unsigned int col = 0; col < state_dimension; ++col) {
                int row_index = row*(number_of_Chebyshev_points-1) + index;
                int col_index = col*(number_of_Chebyshev_points-1) + index;
                // m_C_NN(row_index, col_index) = m_D_NN(row_index, col_index)/upper_integration_limit - A_at_Chebyshev_point(row, col);
                A_NN(row_index, col_index) =  - A_at_Chebyshev_point(row, col);
            }
        }

        index++;
    }

    // std::cout << A_NN << std::endl;




    Eigen::VectorXd diag = A_NN.diagonal();

    const double trace = diag.sum();

    std::cout << "\n\nWith trace : " << trace << "\n\n";
    std::cout.flush();

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp(A_NN);
    auto rank = decomp.rank();

    std::cout << "rank : " << rank << "\n\n";
    std::cout.flush();


    //writeToFile("A_NN", A_NN, "../");

    return 0;
}
