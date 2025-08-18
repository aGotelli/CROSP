#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/numerical_integrators/cosserat_rod_integrators.hpp"
#include "CROSP/numerical_integrators/spectral_method/spectral_integrators.hpp"

#include "utilities/Eigen/eigen_io.hpp"






int main()
{


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
    const unsigned int ne_i = 3;


    ::CROSP::polynomial_representation::PolynomialRepresentation m_polynomial_representation(admitted_deformations, ne_i);


    //  The set of rod properties
    ::CROSP::rod_properties::RodPropertiesSPtr m_rod_properties {
        std::make_shared<::CROSP::rod_properties::RodProperties>(::CROSP::rod_properties::RodDimensions(),
                                                        ::CROSP::rod_properties::MaterialProperties())
    };


    unsigned int m_number_of_Chebyshev_points { 4 };


    ::CROSP::strain_parameterisation_stack::StrainFunction m_constrained_strain {
        ::CROSP::strain_parameterisation_stack::default_constrained_strain
    };


    std::unique_ptr<::CROSP::Spectral> m_cosserat_rod_integrators {
        std::make_unique<::CROSP::Spectral>(m_polynomial_representation,
                                             m_number_of_Chebyshev_points,
                                             m_rod_properties,
                                             m_constrained_strain)
    };


    const unsigned int coordinated_dimension = m_polynomial_representation.getCoordinatesDimension();


    Eigen::VectorXd q = Eigen::VectorXd::Random(coordinated_dimension);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Random(coordinated_dimension);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Random(coordinated_dimension);



    m_cosserat_rod_integrators->updateParameterisation(q, dot_q, ddot_q);




    m_cosserat_rod_integrators->m_idm_integrators->m_quaternion->solveSystem();
    m_cosserat_rod_integrators->m_idm_integrators->m_angular_velocity->solveSystem();


    Eigen::MatrixXd A = m_cosserat_rod_integrators->m_idm_integrators->m_angular_velocity->getC();
    Eigen::MatrixXd b = m_cosserat_rod_integrators->m_idm_integrators->m_angular_velocity->getb();


    const std::string suffix = "_N" + std::to_string(m_number_of_Chebyshev_points) + "m3";
    writeToFile("A" + suffix, A, "tests/data/smw");
    writeToFile("b" + suffix, b, "tests/data/smw");



    //  IMPLEMENTATION FROM GPT WITH SVD DECOMPOSITION

    const unsigned int m = 3;
    const unsigned int n = m_number_of_Chebyshev_points-1;
    const unsigned int N = n*m;

    // Extract the 20x20 diagonal block D from top-left
    Eigen::MatrixXd D = A.block(0, 0, n, n);

    // Construct block-diagonal matrix A0 = blkdiag(D, D, D)
    Eigen::MatrixXd A0 = Eigen::MatrixXd::Zero(N, N);
    for (int i = 0; i < m; ++i) {
        A0.block(i*n, i*n, n, n) = D;
    }

    writeToFile("A0" + suffix, A0, "tests/data/smw");

    // Compute ΔA = A - A0
    Eigen::MatrixXd DeltaA = A - A0;

    writeToFile("DeltaA" + suffix, DeltaA, "tests/data/smw");

    std::cout << "A:\n" << A << "\n\n" << "A0:\n" << A0 << "\n\n" << "DeltaA:\n" << DeltaA << "\n\n";


    Eigen::JacobiSVD<Eigen::MatrixXd> svd(DeltaA, Eigen::ComputeThinU | Eigen::ComputeThinV);

    int r = svd.rank();

    Eigen::VectorXd sigma = svd.singularValues().head(r);
    Eigen::MatrixXd sqrtSigma = sigma.cwiseSqrt().asDiagonal();

    Eigen::MatrixXd U = svd.matrixU().leftCols(r) * sqrtSigma;
    Eigen::MatrixXd V = svd.matrixV().leftCols(r) * sqrtSigma;

    writeToFile("U" + suffix, U, "tests/data/smw");
    writeToFile("V" + suffix, V, "tests/data/smw");
    writeToFile("VT" + suffix, V.transpose(), "tests/data/smw");



//    int blockSize = n;
//    int matrixSize = N;
//    std::vector<Eigen::VectorXd> u_columns;
//    std::vector<Eigen::VectorXd> v_columns;

//    for (int i = 0; i < 3; ++i) {
//        for (int j = 0; j < 3; ++j) {
//            if (i == j) continue;

//            // Extract block B_ij
//            Eigen::MatrixXd Bij = A.block(i * blockSize, j * blockSize, blockSize, blockSize);

//            for (int r = 0; r < blockSize; ++r) {
//                for (int c = 0; c < blockSize; ++c) {
//                    double val = Bij(r, c);
//                    if (val != 0.0) {
//                        Eigen::VectorXd u = Eigen::VectorXd::Zero(matrixSize);
//                        Eigen::VectorXd v = Eigen::VectorXd::Zero(matrixSize);
//                        u(i * blockSize + r) = val;
//                        v(j * blockSize + c) = 1.0;
//                        u_columns.push_back(u);
//                        v_columns.push_back(v);
//                    }
//                }
//            }
//        }
//    }

//    // Build U and V matrices from column vectors
//    int rank = u_columns.size();
//    Eigen::MatrixXd U_v2(matrixSize, rank), V_v2(matrixSize, rank);
//    for (int k = 0; k < rank; ++k) {
//        U.col(k) = u_columns[k];
//        V.col(k) = v_columns[k];
//    }


//    writeToFile("U_v2" + suffix, U_v2, "tests/data/smw");
//    writeToFile("V_v2" + suffix, V_v2, "tests/data/smw");

//    std::cout << "Finished" << std::endl;


    return 0;
}
