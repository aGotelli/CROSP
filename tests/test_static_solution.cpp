#include "CROSP/polynomial_representation/polynomial_representation.hpp"
#include "CROSP/CROSP/cosserat_rod.hpp"

#include "NEMI/newmark_integrator.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include <utilities/Qt/qt_helper.hpp>









int main(int argc, char *argv[])
{
    QApplication* a = new QApplication(argc, argv);

    const std::array<bool, 6> admitted_deformations {
        false,
        true,
        false,
        false,
        false,
        false
    };

    unsigned int n = 4;


    ::CROSP::polynomial_representation::PolynomialRepresentation poly_rep(admitted_deformations);
    ::CROSP::CosseratRod rod(poly_rep);



    Eigen::VectorXd q(::CROSP::polynomial_representation::default_number_of_modes);
    Eigen::VectorXd dot_q, ddot_q;
    dot_q = 0 * q;
    ddot_q = 0 * q;

    rod.updateParameterisation(q, dot_q, ddot_q);



    ::LieAlgebra::Vector6d F1;
    F1   << 0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0;

    q.setZero();
    CROSP::static_initialiser::Initializer initialiser(rod, F1, q);


    ::Chebyshev::ChebyshevReconstructor Chebyshev_reconstructor;


    const auto stack_of_points = Chebyshev_reconstructor.reconstructRodShape(
                                        rod.getRodPositionsAtChebyshevPoints());

    ::utilities::qt_helper::plot2D(a,
                                   stack_of_points(Eigen::all, {0, 2}),
                                   "Rod Static Solution",
                                   "x [m]",
                                   "z [m]",
                                   "Rod Shape");



    return 0;
}
