#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/idm_integrators/idm_integrators.hpp"
#include "CROSP/CROSP/cosserat_rod.hpp"
#include "CROSP/tidm_integrators/tidm_integrators.hpp"
#include "CROSP/polynomial_representation/polynomial_representation.hpp"


#include <utilities/Eigen/eigen_io.hpp>

#include <benchmark/benchmark.h>

#include <fstream>

using namespace CROSP;




int main(int argc, char *argv[])
{

    //  =========================================================================================
    //                               PARAMETERS FOR THE INTEGRATIONS
    //  =========================================================================================

    const unsigned int number_of_Chebyshev_points = 17;

    //  Number of modes per deformation
    const unsigned int ne = 3;

    //  Deformation you want to allow
    const std::array<bool, 6> admitted_deformations {
        false,  //  torsion along x
        true,   //  bendin along y
        false,  //  benging along z
        false,  //  elongation along x
        false,  //  shear on y
        false   //  shear on z
    };

    const double radius = 0.001; //mm
    const double length = 1.0; // m


    const std::string path_where_to_write_the_results = "data/"; //  If empty it will write in the build folder
    //  Note that if it starts with a / then the path is absolute (starting from the home directory of your pc)
    //  Otherwise the path is relative to your build directory.
    //  Note that the directory you want to save the data in should exist before running this code


    auto polynomial_representation =
            std::make_shared<::polynomial_representation::PolynomialRepresentation>(admitted_deformations,
                                                                                   ne);

    Eigen::VectorXd qe = Eigen::VectorXd::Zero(polynomial_representation->getCoordinatesDimension());

    //  Generalised coordinates
    qe << -0.2,
           0.2,
          -0.2;

    //  Define the wrench to apply at the tip (X=1) in the coordinates of the reference frame (world coordinates)
    const Eigen::Vector3d N1(0, 0.5, -1.5);
    const Eigen::Vector3d C1(0, 0, -1.0);



    //  =========================================================================================
    //                                  PREPROCESSING
    //  =========================================================================================

    ::CROSP::rod_properties::CircularCrossSection cs(radius);
    ::CROSP::rod_properties::RodDimensions rod_dimensions(&cs, length);


    auto strain_param =
            std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(polynomial_representation,
                                                                                       number_of_Chebyshev_points);


    std::shared_ptr<::rod_properties::RodProperties> rod_properties =
            std::make_shared<::rod_properties::RodProperties>(strain_param,
                                                              rod_dimensions,
                                                              ::CROSP::rod_properties::MaterialProperties());



    const std::shared_ptr<::strain_parameterisation::StrainParameterisation> strain_parameterisation =
            std::make_shared<::strain_parameterisation::StrainParameterisation>(polynomial_representation,
                                                                                number_of_Chebyshev_points);

    std::shared_ptr<::strain_parameterisation::StrainParameterisation> strain_parameterisation_Delta =
            std::make_shared<::strain_parameterisation::StrainParameterisation>(strain_parameterisation,
                                                                                ::LieAlgebra::Vector6d::Zero());



    std::shared_ptr<::idm_integrators::IDMIntegrators> idm_integrators =
            std::make_shared<::idm_integrators::IDMIntegrators>(strain_parameterisation,
                                                                rod_properties);






    Eigen::VectorXd dot_qe = Eigen::VectorXd::Zero(polynomial_representation->getCoordinatesDimension());
    Eigen::VectorXd ddot_qe = Eigen::VectorXd::Zero(polynomial_representation->getCoordinatesDimension());


    //  =========================================================================================
    //                                  COMPUTATIONS
    //  =========================================================================================

    //  Update strains
    strain_parameterisation->updateStacks(qe, dot_qe, ddot_qe);


    //  All the integrations for the forward
    idm_integrators->m_quaternion->integrate( Eigen::Vector4d(1, 0, 0, 0) );

    idm_integrators->m_position->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_linear_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_acceleration->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_linear_acceleration->integrate( Eigen::Vector3d::Zero() );

    //  Map force and couple into local coordinates
    const auto q1 = idm_integrators->m_quaternion->getStateAtPoint(0);
    const auto R1 = Eigen::Quaterniond(q1[0], q1[1], q1[2], q1[3]).toRotationMatrix();
    Eigen::Vector3d force_at_tip_local_coord = R1.transpose()*N1;
    Eigen::Vector3d couple_at_tip_local_coord = R1.transpose()*C1;


    //  integrations for the bckward
    idm_integrators->m_internal_forces->integrate( force_at_tip_local_coord );

    idm_integrators->m_internal_couples->integrate( couple_at_tip_local_coord );

    idm_integrators->m_generalised_forces->integrate( Eigen::VectorXd::Zero(polynomial_representation->m_total_number_of_modes) );

    std::cout.flush();


    //  =========================================================================================
    //                                  SAVE DATA
    //  =========================================================================================

    writeToFile("Quaternions", idm_integrators->m_quaternion->getStackAsMatrix(), path_where_to_write_the_results);
    writeToFile("Positions", idm_integrators->m_position->getStackAsMatrix(), path_where_to_write_the_results);
    writeToFile("AngularVelocities", idm_integrators->m_angular_velocity->getStackAsMatrix(), path_where_to_write_the_results);
    writeToFile("LinearVelocities", idm_integrators->m_linear_velocity->getStackAsMatrix(), path_where_to_write_the_results);
    writeToFile("AngularAccelerations", idm_integrators->m_angular_acceleration->getStackAsMatrix(), path_where_to_write_the_results);
    writeToFile("LinearAccelerations", idm_integrators->m_linear_acceleration->getStackAsMatrix(), path_where_to_write_the_results);
    writeToFile("InternalCouples", idm_integrators->m_internal_couples->getStackAsMatrix(), path_where_to_write_the_results);
    writeToFile("InternalForces", idm_integrators->m_internal_forces->getStackAsMatrix(), path_where_to_write_the_results);
    writeToFile("GeneralisedForces", idm_integrators->m_generalised_forces->getStackAsMatrix(), path_where_to_write_the_results);


    //  =========================================================================================
    //                                  LOG DATA
    //  =========================================================================================

    std::cout << "Quaternions : \n" << idm_integrators->m_quaternion->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Positions : \n" << idm_integrators->m_position->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Angular velocities : \n" << idm_integrators->m_angular_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Linear velocities : \n" << idm_integrators->m_linear_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Angular accelerations : \n" << idm_integrators->m_angular_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Linear accelerations : \n" << idm_integrators->m_linear_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;

    std::cout << "Internal couples : \n" << idm_integrators->m_internal_couples->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Internal forces : \n" << idm_integrators->m_internal_forces->getStackAsMatrix() << "\n\n\n" << std::endl;






    return 0;
}
