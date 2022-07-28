#include "strain_parameterisation.hpp"
#include "idm_integrators.hpp"
#include "cosserat_rod.hpp"
#include "tidm_integrators.hpp"



using namespace CROSP;




int main(int argc, char *argv[])
{

    const unsigned int number_of_Chebyshev_points = 17;

    const unsigned int ne = 3;
    const unsigned int na = 1;


    std::shared_ptr<StrainParameterisation> strain_parameterisation = std::make_shared<StrainParameterisation>(ne, na,
                                                                                                               number_of_Chebyshev_points);

    std::shared_ptr<StrainParameterisationPerturbation> strain_parameterisation_perturbation = std::make_shared<StrainParameterisationPerturbation>(ne, na,
                                                                                                                                                    number_of_Chebyshev_points);

    Eigen::Matrix3d M_angular, M_linear;
    M_angular.setIdentity();
    M_linear.setIdentity();

    std::shared_ptr<IDMIntegrators> idm_integrators = std::make_shared<IDMIntegrators>(strain_parameterisation, M_angular, M_linear);

    std::shared_ptr<TIDMIntegrators> tidm_integrators = std::make_shared<TIDMIntegrators>(strain_parameterisation,
                                                                                          strain_parameterisation_perturbation,
                                                                                          idm_integrators);

    Eigen::VectorXd qe = Eigen::VectorXd::Zero(ne*na);
    Eigen::VectorXd dot_qe = Eigen::VectorXd::Zero(ne*na);
    Eigen::VectorXd ddot_qe = Eigen::VectorXd::Zero(ne*na);

    strain_parameterisation->update(qe, dot_qe, ddot_qe);

    idm_integrators->m_quaternion->integrate( Eigen::Vector4d(1, 0, 0, 0) );

    idm_integrators->m_position->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_linear_velocity->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_angular_acceleration->integrate( Eigen::Vector3d::Zero() );

    idm_integrators->m_linear_acceleration->integrate( Eigen::Vector3d::Zero() );


    Eigen::VectorXd Delta_qe = Eigen::VectorXd::Zero(ne*na);
    Eigen::VectorXd Delta_dot_qe = Eigen::VectorXd::Zero(ne*na);
    Eigen::VectorXd Delta_ddot_qe = Eigen::VectorXd::Zero(ne*na);


    const unsigned int a = 400;
    const unsigned int b = 160000;

    for(unsigned int i=0; i<3; i++){
        Delta_qe.setZero();
        Delta_qe(i) = 1;

        Delta_dot_qe = a*Delta_qe;
        Delta_ddot_qe = b*Delta_qe;

        strain_parameterisation_perturbation->update(Delta_qe, Delta_dot_qe, Delta_ddot_qe);

        tidm_integrators->m_delta_rotation->integrate( Eigen::Vector3d::Zero() );
        tidm_integrators->m_delta_position->integrate( Eigen::Vector3d::Zero() );


        Eigen::Matrix<double, 6, number_of_Chebyshev_points> Delta_zeta;
        for(unsigned int point=0; point<number_of_Chebyshev_points; point++){
            Delta_zeta.block<3,1>(0, point) = tidm_integrators->m_delta_rotation->getStateAtPoint(point);
            Delta_zeta.block<3,1>(3, point) = tidm_integrators->m_delta_position->getStateAtPoint(point);
        }

        std::cout << Delta_zeta << std::endl;
    }





    return 0;
}
