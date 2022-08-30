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

    const unsigned int number_of_Chebyshev_points = 21;


    const unsigned int ne = 3;
    const std::array<bool, 6> admitted_deformations {
        false,
        true,
        false,
        false,
        false,
        false
    };


    ::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations,
                                                                                    ne);


    ::CROSP::rod_properties::RodDimensions rod_dimensions(0.01, 1.0);
    ::CROSP::rod_properties::MaterialProperties material_properties(210e9, 80e9, 0);

    std::shared_ptr<::rod_properties::RodProperties> rod_properties =
            std::make_shared<::rod_properties::RodProperties>(polynomial_representation,
//                                                              material_properties,
                                                              rod_dimensions);



    std::shared_ptr<::strain_parameterisation::StrainParameterisation> strain_parameterisation =
            std::make_shared<::strain_parameterisation::StrainParameterisation>(polynomial_representation,
                                                                                number_of_Chebyshev_points);

    std::shared_ptr<::strain_parameterisation::StrainParameterisation> strain_parameterisation_Delta =
            std::make_shared<::strain_parameterisation::StrainParameterisation>(polynomial_representation,
                                                                                ::LieAlgebra::Vector6d::Zero(),
                                                                                number_of_Chebyshev_points);



    std::shared_ptr<::idm_integrators::IDMIntegrators> idm_integrators =
            std::make_shared<::idm_integrators::IDMIntegrators>(number_of_Chebyshev_points,
                                                                polynomial_representation,
                                                                strain_parameterisation,
                                                                rod_properties);

    std::shared_ptr<::tidm_integrators::TIDMIntegrators> tidm_integrators =
            std::make_shared<::tidm_integrators::TIDMIntegrators>(number_of_Chebyshev_points,
                                                                  polynomial_representation,
                                                                  strain_parameterisation,
                                                                  strain_parameterisation_Delta,
                                                                  idm_integrators,
                                                                  rod_properties);

    Eigen::VectorXd qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());
    qe << -0.2,
           0.2,
          -0.2;

    Eigen::VectorXd dot_qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());
    dot_qe = 0.4 * qe;
    Eigen::VectorXd ddot_qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());
    ddot_qe = -0.15*qe;


    const Eigen::Vector3d N1(0, 0.5, -1.5);
    const Eigen::Vector3d C1(0, 0, -1.0);

    strain_parameterisation->updateStacks(qe, dot_qe, ddot_qe);

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

    idm_integrators->m_internal_forces->integrate( force_at_tip_local_coord );

    idm_integrators->m_internal_couples->integrate( couple_at_tip_local_coord );

    idm_integrators->m_generalised_forces->integrate( Eigen::VectorXd::Zero(polynomial_representation.m_ne) );

    std::cout.flush();

    writeToFile("Quaternions", idm_integrators->m_quaternion->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
    writeToFile("Positions", idm_integrators->m_position->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
    writeToFile("AngularVelocities", idm_integrators->m_angular_velocity->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
    writeToFile("LinearVelocities", idm_integrators->m_linear_velocity->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
    writeToFile("AngularAccelerations", idm_integrators->m_angular_acceleration->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
    writeToFile("LinearAccelerations", idm_integrators->m_linear_acceleration->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
    writeToFile("InternalCouples", idm_integrators->m_internal_couples->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
    writeToFile("InternalForces", idm_integrators->m_internal_forces->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
    writeToFile("GeneralisedForces", idm_integrators->m_generalised_forces->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");



    std::cout << "Quaternions : \n" << idm_integrators->m_quaternion->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Positions : \n" << idm_integrators->m_position->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Angular velocities : \n" << idm_integrators->m_angular_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Linear velocities : \n" << idm_integrators->m_linear_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Angular accelerations : \n" << idm_integrators->m_angular_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Linear accelerations : \n" << idm_integrators->m_linear_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;

    std::cout << "Internal couples : \n" << idm_integrators->m_internal_couples->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Internal forces : \n" << idm_integrators->m_internal_forces->getStackAsMatrix() << "\n\n\n" << std::endl;


    Eigen::VectorXd Delta_qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());
    Eigen::VectorXd Delta_dot_qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());
    Eigen::VectorXd Delta_ddot_qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());


    const unsigned int a = 400;
    const unsigned int b = 160000;




    Eigen::MatrixXd Delta_zeta(polynomial_representation.m_ne*6, number_of_Chebyshev_points);
    for(unsigned int i=0; i<polynomial_representation.m_ne; i++){
        Delta_qe.setZero();
        Delta_qe(i) = 1;

        Delta_dot_qe = a*Delta_qe;
        Delta_ddot_qe = b*Delta_qe;

        strain_parameterisation_Delta->updateStacks(Delta_qe, Delta_dot_qe, Delta_ddot_qe);

        tidm_integrators->m_Delta_rotation->integrate( Eigen::Vector3d::Zero() );
        tidm_integrators->m_Delta_position->integrate( Eigen::Vector3d::Zero() );



        std::cout << "Delta : " << (i+1) << "\n\n";
//        std::cout << "  - Delta rotations : \n" << tidm_integrators->m_Delta_rotation->getStackAsMatrix() << "\n";
//        std::cout << "  - Delta positions : \n" << tidm_integrators->m_Delta_position->getStackAsMatrix() << "\n";
        writeToFile("Delta_P"+std::to_string(i+1), tidm_integrators->m_Delta_rotation->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
        writeToFile("Delta_Y"+std::to_string(i+1), tidm_integrators->m_Delta_position->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");


//        std::cout << "\n\n\n";

        tidm_integrators->m_Delta_angular_velocity->integrate(Eigen::Vector3d::Zero());
        tidm_integrators->m_Delta_linear_velocity->integrate(Eigen::Vector3d::Zero());
//        std::cout << "  - Delta angular velocity : \n" << tidm_integrators->m_Delta_angular_velocity->getStackAsMatrix() << "\n";
//        std::cout << "  - Delta linear velocity : \n" << tidm_integrators->m_Delta_linear_velocity->getStackAsMatrix() << "\n";
        writeToFile("Delta_Omega"+std::to_string(i+1), tidm_integrators->m_Delta_angular_velocity->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
        writeToFile("Delta_V"+std::to_string(i+1), tidm_integrators->m_Delta_linear_velocity->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");

//        std::cout << "\n\n\n";


        tidm_integrators->m_Delta_angular_acceleration->integrate(Eigen::Vector3d::Zero());
        tidm_integrators->m_Delta_linear_acceleration->integrate(Eigen::Vector3d::Zero());
        writeToFile("Delta_dot_Omega"+std::to_string(i+1), tidm_integrators->m_Delta_angular_acceleration->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
        writeToFile("Delta_dot_V"+std::to_string(i+1), tidm_integrators->m_Delta_linear_acceleration->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");



//        std::cout << "\n\n\n";


        tidm_integrators->m_Delta_internal_forces->integrate(Eigen::Vector3d::Zero());
        tidm_integrators->m_Delta_internal_couples->integrate(Eigen::Vector3d::Zero());
        std::cout << "  - Delta internal forces : \n" << tidm_integrators->m_Delta_internal_forces->getStackAsMatrix() << "\n";
        std::cout << "  - Delta internal couples : \n" << tidm_integrators->m_Delta_internal_couples->getStackAsMatrix() << "\n";

        writeToFile("Delta_N"+std::to_string(i+1), tidm_integrators->m_Delta_internal_forces->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");
        writeToFile("Delta_C"+std::to_string(i+1), tidm_integrators->m_Delta_internal_couples->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");


        tidm_integrators->m_Delta_generalised_forces->integrate(Eigen::VectorXd::Zero(polynomial_representation.m_ne));
        writeToFile("Delta_Q"+std::to_string(i+1), tidm_integrators->m_Delta_generalised_forces->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");

        std::cout << "\n\n\n\n\n\n";

    }

    ::benchmark::RegisterBenchmark("Compose Jacobian", [&](::benchmark::State &t_state){

        while(t_state.KeepRunning()){

            strain_parameterisation->updateStacks(qe, dot_qe, ddot_qe);

            idm_integrators->m_quaternion->integrate( Eigen::Vector4d(1, 0, 0, 0) );

            idm_integrators->m_position->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_angular_velocity->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_linear_velocity->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_angular_acceleration->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_linear_acceleration->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_internal_forces->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_internal_couples->integrate( Eigen::Vector3d::Zero() );


            idm_integrators->m_generalised_forces->integrate( Eigen::VectorXd::Zero(polynomial_representation.m_ne) );


            for(unsigned int i=0; i<polynomial_representation.m_ne; i++){
                Delta_qe.setZero();
                Delta_qe(i) = 1;

                Delta_dot_qe = a*Delta_qe;
                Delta_ddot_qe = b*Delta_qe;

                strain_parameterisation_Delta->updateStacks(Delta_qe, Delta_dot_qe, Delta_ddot_qe);

                tidm_integrators->m_Delta_rotation->integrate( Eigen::Vector3d::Zero() );
                tidm_integrators->m_Delta_position->integrate( Eigen::Vector3d::Zero() );


                tidm_integrators->m_Delta_angular_velocity->integrate(Eigen::Vector3d::Zero());
                tidm_integrators->m_Delta_linear_velocity->integrate(Eigen::Vector3d::Zero());


                tidm_integrators->m_Delta_angular_acceleration->integrate(Eigen::Vector3d::Zero());
                tidm_integrators->m_Delta_linear_acceleration->integrate(Eigen::Vector3d::Zero());

                tidm_integrators->m_Delta_internal_forces->integrate(Eigen::Vector3d::Zero());
                tidm_integrators->m_Delta_internal_couples->integrate(Eigen::Vector3d::Zero());

            }

        }
    })->Repetitions(50);

    ::benchmark::Initialize(&argc, argv);

    //::benchmark::RunSpecifiedBenchmarks();





    return 0;
}
