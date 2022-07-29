#include "strain_parameterisation.hpp"
#include "idm_integrators.hpp"
#include "cosserat_rod.hpp"
#include "tidm_integrators.hpp"



using namespace CROSP;




int main(int argc, char *argv[])
{

    const unsigned int number_of_Chebyshev_points = 11;

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



    std::shared_ptr<StrainParameterisation> strain_parameterisation = std::make_shared<StrainParameterisation>(ne, admitted_deformations, number_of_Chebyshev_points);

    std::shared_ptr<StrainParameterisationPerturbation> strain_parameterisation_perturbation = std::make_shared<StrainParameterisationPerturbation>(ne, admitted_deformations,
                                                                                                                                                    number_of_Chebyshev_points);

    Eigen::Matrix3d M_angular, M_linear;
    M_angular.setIdentity();
    M_linear.setIdentity();

    std::shared_ptr<IDMIntegrators> idm_integrators = std::make_shared<IDMIntegrators>(strain_parameterisation,
                                                                                       M_angular, M_linear,
                                                                                       number_of_Chebyshev_points);

    std::shared_ptr<TIDMIntegrators> tidm_integrators = std::make_shared<TIDMIntegrators>(strain_parameterisation,
                                                                                          strain_parameterisation_perturbation,
                                                                                          idm_integrators,
                                                                                          number_of_Chebyshev_points);

    Eigen::VectorXd qe = Eigen::VectorXd::Zero(ne*na);
    qe << -0.2,
           0.2,
          -0.2;

    Eigen::VectorXd dot_qe = Eigen::VectorXd::Zero(ne*na);
    //dot_qe = 0.4 * qe;
    Eigen::VectorXd ddot_qe = Eigen::VectorXd::Zero(ne*na);
    //ddot_qe = -0.15*qe;

    strain_parameterisation->update(qe, dot_qe, ddot_qe);

    idm_integrators->m_quaternion->integrate( Eigen::Vector4d(1, 0, 0, 0) );

    idm_integrators->m_position->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_linear_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_acceleration->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_linear_acceleration->integrate( Eigen::Vector3d::Zero() );

//    std::cout << "Quaternions : \n" << idm_integrators->m_quaternion->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Positions : \n" << idm_integrators->m_position->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Angular velocities : \n" << idm_integrators->m_angular_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Linear velocities : \n" << idm_integrators->m_linear_velocity->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Angular accelerations : \n" << idm_integrators->m_angular_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;
//    std::cout << "Linear accelerations : \n" << idm_integrators->m_linear_acceleration->getStackAsMatrix() << "\n\n\n" << std::endl;


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
        std::cout << "  - rotations : \n" << tidm_integrators->m_Delta_rotation->getStackAsMatrix() << "\n";
        std::cout << "  - positions : \n" << tidm_integrators->m_Delta_position->getStackAsMatrix() << "\n";

        std::cout << "\n\n\n";

        tidm_integrators->m_Delta_angular_velocity->integrate(Eigen::Vector3d::Zero());
        tidm_integrators->m_Delta_linear_velocity->integrate(Eigen::Vector3d::Zero());
        std::cout << "  - angular velocity : \n" << tidm_integrators->m_Delta_angular_velocity->getStackAsMatrix() << "\n";
        std::cout << "  - linear velocity : \n" << tidm_integrators->m_Delta_linear_velocity->getStackAsMatrix() << "\n";

        std::cout << "\n\n\n";


        tidm_integrators->m_Delta_angular_acceleration->integrate(Eigen::Vector3d::Zero());
//        tidm_integrators->m_Delta_linear_acceleration->integrate(Eigen::Vector3d::Zero());
        std::cout << "  - angular acceleration : \n" << tidm_integrators->m_Delta_angular_acceleration->getStackAsMatrix() << "\n";
//        std::cout << "  - linear acceleration : \n" << tidm_integrators->m_Delta_linear_acceleration->getStackAsMatrix() << "\n";


//        Eigen::Matrix<double, 6, number_of_Chebyshev_points> local_Delta_zeta;
//        for(unsigned int point=0; point<number_of_Chebyshev_points; point++){




//            local_Delta_zeta.block<3,1>(0, point) = tidm_integrators->m_delta_rotation->getStateAtPoint(point);
//            local_Delta_zeta.block<3,1>(3, point) = tidm_integrators->m_delta_position->getStateAtPoint(point);

//            std::cout << "At delta : " << point << " Delta zeta \n" << ( Eigen::VectorXd(6) << tidm_integrators->m_delta_rotation->getStateAtPoint(point), tidm_integrators->m_delta_position->getStateAtPoint(point) ).finished() << "\n\n";

//            //auto row = i*6;
//            //Delta_zeta.block<6, number_of_Chebyshev_points>(row + 0, point) = local_Delta_zeta;
//        }
//        std::cout << "At delta : " << i << " Delta zeta \n" << local_Delta_zeta.reverse() << "\n\n";

    }

    //std::cout << "All the Delta zeta : \n\n" << Delta_zeta << std::endl;





    return 0;
}
