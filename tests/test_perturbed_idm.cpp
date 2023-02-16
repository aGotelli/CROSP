#include "CROSP/strain_parameterisation/strain_parameterisation.hpp"
#include "CROSP/idm_integrators/idm_integrators.hpp"
#include "CROSP/CROSP/cosserat_rod.hpp"
#include "CROSP/tidm_integrators/tidm_integrators.hpp"
#include "CROSP/polynomial_representation/polynomial_representation.hpp"


#include <utilities/Eigen/eigen_io.hpp>


#include <fstream>

using namespace CROSP;




int main(int argc, char *argv[])
{

    const std::array<bool, 6> admitted_deformations {
        false,
        true,
        false,
        false,
        false,
        false
    };

    const unsigned int ne = 4;

    const unsigned int na = std::count(admitted_deformations.begin(),
                                                  admitted_deformations.end(),
                                                  true);

    const unsigned int coordinates_dimension = ne * na;

    const unsigned int number_of_Chebyshev_points = 17;





    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    const double mu = 0.0;
    std::shared_ptr<::rod_properties::RodProperties> rod_properties =
            std::make_shared<::rod_properties::RodProperties>(mu,
                                                              polynomial_representation);

    std::shared_ptr<::strain_parameterisation::StrainParameterisation> strain_parameterisation =
            std::make_shared<::strain_parameterisation::StrainParameterisation>(polynomial_representation,
                                                                                number_of_Chebyshev_points);

    std::shared_ptr<::idm_integrators::IDMIntegrators> idm_integrators =
            std::make_shared<::idm_integrators::IDMIntegrators>(number_of_Chebyshev_points,
                                                                polynomial_representation,
                                                                strain_parameterisation,
                                                                rod_properties);





    Eigen::VectorXd qe = Eigen::VectorXd::Zero(coordinates_dimension);
    Eigen::VectorXd dot_qe = Eigen::VectorXd::Zero(coordinates_dimension);
    Eigen::VectorXd ddot_qe = Eigen::VectorXd::Zero(coordinates_dimension);


    const Eigen::Vector3d N1 = Eigen::Vector3d::Zero();
    const Eigen::Vector3d C1 = Eigen::Vector3d::Zero();

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

    idm_integrators->m_generalised_forces->integrate( Eigen::VectorXd::Zero(coordinates_dimension) );


    std::cout << "Quaternions : \n" << idm_integrators->m_quaternion->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Positions : \n" << idm_integrators->m_position->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Angular velocities : \n" << idm_integrators->m_angular_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Linear velocities : \n" << idm_integrators->m_linear_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Angular accelerations : \n" << idm_integrators->m_angular_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Linear accelerations : \n" << idm_integrators->m_linear_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;

    std::cout << "Internal couples : \n" << idm_integrators->m_internal_couples->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Internal forces : \n" << idm_integrators->m_internal_forces->getStackAsMatrix() << "\n\n\n" << std::endl;
    std::cout << "Generalised forces : \n" << idm_integrators->m_generalised_forces->getStackAsMatrix() << "\n\n\n" << std::endl;

    Eigen::VectorXd Qe = idm_integrators->m_generalised_forces->getStateAtPoint(OSNI::INTEGRATION_DOMAIN::BEGIN);
    std::cout << "Origianl Qe : \n" << Qe << "\n\n\n" << std::endl;

    Eigen::VectorXd elastic_energy = rod_properties->m_Kee * qe;
    std::cout << "Origianl elastic energy : \n" << elastic_energy << "\n\n\n" << std::endl;

    Eigen::VectorXd residual = elastic_energy - Qe;
    std::cout << "Origianl residual : \n" << residual << "\n\n\n" << std::endl;

    std::cout.flush();






    const unsigned int a = 400;
    const unsigned int b = 160000;


    Eigen::VectorXd delta_qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());
    Eigen::VectorXd delta_dot_qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());
    Eigen::VectorXd delta_ddot_qe = Eigen::VectorXd::Zero(polynomial_representation.getCoordinatesDimension());

    Eigen::VectorXd q_variation = Eigen::VectorXd::Zero(coordinates_dimension);
    const double delta = 1e-6;

    for(unsigned int i=0; i<polynomial_representation.m_total_number_of_modes; i++){

        std::cout << "\n\n\n\n\n\n\n\n\n" << "i : " << i << "\n\n\n\n\n\n";

        delta_qe = qe;
        delta_qe(i) += delta;


        delta_dot_qe = a*delta_qe;
        delta_ddot_qe = b*delta_qe;

        strain_parameterisation->updateStacks(delta_qe, delta_dot_qe, delta_ddot_qe);

        idm_integrators->m_quaternion->integrate( Eigen::Vector4d(1, 0, 0, 0) );
        idm_integrators->m_position->integrate( Eigen::Vector3d::Zero() );

        idm_integrators->m_angular_velocity->integrate( Eigen::Vector3d::Zero() );
        idm_integrators->m_linear_velocity->integrate( Eigen::Vector3d::Zero() );

        idm_integrators->m_angular_acceleration->integrate( Eigen::Vector3d::Zero() );
        idm_integrators->m_linear_acceleration->integrate( Eigen::Vector3d::Zero() );

        //  Map force and couple into local coordinates
        const auto delta_q1 = idm_integrators->m_quaternion->getStateAtPoint(0);
        const auto delta_R1 = Eigen::Quaterniond(delta_q1[0],
                                                 delta_q1[1],
                                                 delta_q1[2],
                                                 delta_q1[3]).toRotationMatrix();
        Eigen::Vector3d force_at_tip_local_coord = delta_R1.transpose()*N1;
        Eigen::Vector3d couple_at_tip_local_coord = delta_R1.transpose()*C1;

        idm_integrators->m_internal_forces->integrate( force_at_tip_local_coord );
        idm_integrators->m_internal_couples->integrate( couple_at_tip_local_coord );

        idm_integrators->m_generalised_forces->integrate( Eigen::VectorXd::Zero(coordinates_dimension) );


        std::cout << "Quaternions : \n" << idm_integrators->m_quaternion->getStackAsMatrix() << "\n\n\n" << std::endl;
        std::cout << "Positions : \n" << idm_integrators->m_position->getStackAsMatrix() << "\n\n\n" << std::endl;
        std::cout << "Angular velocities : \n" << idm_integrators->m_angular_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
        std::cout << "Linear velocities : \n" << idm_integrators->m_linear_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
        std::cout << "Angular accelerations : \n" << idm_integrators->m_angular_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;
        std::cout << "Linear accelerations : \n" << idm_integrators->m_linear_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;

        std::cout << "Internal couples : \n" << idm_integrators->m_internal_couples->getStackAsMatrix() << "\n\n\n" << std::endl;
        std::cout << "Internal forces : \n" << idm_integrators->m_internal_forces->getStackAsMatrix() << "\n\n\n" << std::endl;
        std::cout << "Generalised forces : \n" << idm_integrators->m_generalised_forces->getStackAsMatrix() << "\n\n\n" << std::endl;


        std::cout << "Origianl Qe : \n" << Qe << "\n" << std::endl;

        Eigen::VectorXd delta_Qe = idm_integrators->m_generalised_forces->getStateAtPoint(OSNI::INTEGRATION_DOMAIN::BEGIN);

        Eigen::VectorXd delta_elastic_energy = rod_properties->m_Kee * delta_qe;
        std::cout << "delta elastic energy : \n" << delta_elastic_energy << "\n" << std::endl;

        Eigen::VectorXd delta_residual = delta_elastic_energy - delta_Qe;
        std::cout << "delta residual : \n" << delta_residual << "\n" << std::endl;

        std::cout << "Origianl residual : \n" << residual << "\n" << std::endl;

        Eigen::VectorXd residual_variation = delta_residual - residual;
        std::cout << "delta_residual - residual : \n" << residual_variation << "\n" << std::endl;

        Eigen::VectorXd col = residual_variation / delta;
        std::cout << "Jacobian column : \n" << col << "\n" << std::endl;


        std::cout.flush();

    }





    return 0;
}
