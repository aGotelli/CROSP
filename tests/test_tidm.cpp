#include "strain_parameterisation.hpp"
#include "idm_integrators.hpp"
#include "cosserat_rod.hpp"
#include "tidm_integrators.hpp"

#include <utilities/Eigen/eigen_io.hpp>

#include <benchmark/benchmark.h>

#include <fstream>

using namespace CROSP;





int main(int argc, char *argv[])
{

    const unsigned int number_of_Chebyshev_points = 21;

    const unsigned int ne = 3;


    const unsigned int na = 1;
    std::vector<bool> admitted_deformations {
        false,
        true,
        false,
        false,
        false,
        false
    };



    std::shared_ptr<MaterialProperties> material_properties = std::make_shared<MaterialProperties>();



    std::shared_ptr<StrainParameterisation> strain_parameterisation = std::make_shared<StrainParameterisation>(ne,
                                                                                                               admitted_deformations,
                                                                                                               number_of_Chebyshev_points);

    std::shared_ptr<StrainParameterisationPerturbation> strain_parameterisation_perturbation = std::make_shared<StrainParameterisationPerturbation>(ne,
                                                                                                                                                    admitted_deformations,
                                                                                                                                                    number_of_Chebyshev_points);



    std::shared_ptr<IDMIntegrators> idm_integrators = std::make_shared<IDMIntegrators>(strain_parameterisation,
                                                                                       material_properties,
                                                                                       number_of_Chebyshev_points);

    std::shared_ptr<TIDMIntegrators> tidm_integrators = std::make_shared<TIDMIntegrators>(strain_parameterisation,
                                                                                          strain_parameterisation_perturbation,
                                                                                          idm_integrators,
                                                                                          material_properties,
                                                                                          number_of_Chebyshev_points);

    Eigen::VectorXd qe = Eigen::VectorXd::Zero(ne*na);
    qe << -0.2,
           0.2,
          -0.2;

    Eigen::VectorXd dot_qe = Eigen::VectorXd::Zero(ne*na);
    dot_qe = 0.4 * qe;
    Eigen::VectorXd ddot_qe = Eigen::VectorXd::Zero(ne*na);
    ddot_qe = -0.15*qe;

    strain_parameterisation->update(qe, dot_qe, ddot_qe);

    idm_integrators->m_quaternion->integrate( Eigen::Vector4d(1, 0, 0, 0) );

    idm_integrators->m_position->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_linear_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_acceleration->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_linear_acceleration->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_internal_forces->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_internal_couples->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_generalised_forces->integrate( Eigen::VectorXd::Zero(ne) );

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



//    std::cout << "Quaternions : \n" << idm_integrators->m_quaternion->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Positions : \n" << idm_integrators->m_position->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Angular velocities : \n" << idm_integrators->m_angular_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Linear velocities : \n" << idm_integrators->m_linear_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Angular accelerations : \n" << idm_integrators->m_angular_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Linear accelerations : \n" << idm_integrators->m_linear_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;

//    std::cout << "Internal couples : \n" << idm_integrators->m_internal_couples->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Internal forces : \n" << idm_integrators->m_internal_forces->getStackAsMatrix() << "\n\n\n" << std::endl;


    Eigen::VectorXd Delta_qe = Eigen::VectorXd::Zero(ne*na);
    Eigen::VectorXd Delta_dot_qe = Eigen::VectorXd::Zero(ne*na);
    Eigen::VectorXd Delta_ddot_qe = Eigen::VectorXd::Zero(ne*na);


    const unsigned int a = 400;
    const unsigned int b = 160000;




    Eigen::Matrix<double, ne*6, number_of_Chebyshev_points> Delta_zeta;
    for(unsigned int i=0; i<ne; i++){
        Delta_qe.setZero();
        Delta_qe(i) = 1;

        Delta_dot_qe = a*Delta_qe;
        Delta_ddot_qe = b*Delta_qe;

        strain_parameterisation_perturbation->update(Delta_qe, Delta_dot_qe, Delta_ddot_qe);

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


        tidm_integrators->m_Delta_generalised_forces->integrate(Eigen::VectorXd::Zero(ne));
        writeToFile("Delta_Q"+std::to_string(i+1), tidm_integrators->m_Delta_generalised_forces->getStackAsMatrix(), "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");

        std::cout << "\n\n\n\n\n\n";

    }

    ::benchmark::RegisterBenchmark("Compose Jacobian", [&](::benchmark::State &t_state){

        while(t_state.KeepRunning()){

            strain_parameterisation->update(qe, dot_qe, ddot_qe);

            idm_integrators->m_quaternion->integrate( Eigen::Vector4d(1, 0, 0, 0) );

            idm_integrators->m_position->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_angular_velocity->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_linear_velocity->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_angular_acceleration->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_linear_acceleration->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_internal_forces->integrate( Eigen::Vector3d::Zero() );

            idm_integrators->m_internal_couples->integrate( Eigen::Vector3d::Zero() );


            idm_integrators->m_generalised_forces->integrate( Eigen::VectorXd::Zero(ne) );


            for(unsigned int i=0; i<ne; i++){
                Delta_qe.setZero();
                Delta_qe(i) = 1;

                Delta_dot_qe = a*Delta_qe;
                Delta_ddot_qe = b*Delta_qe;

                strain_parameterisation_perturbation->update(Delta_qe, Delta_dot_qe, Delta_ddot_qe);

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
