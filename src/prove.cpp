#include "CROSP/CROSP/cosserat_rod.hpp"




#include "ATORS/simple_actuator/simple_actuator.hpp"
#include "ATORS/tendon_driven_actuation/tendon_driven_actuation.hpp"
#include "CROSP/CROSP/internally_actuated_cosserat_rod.hpp"








int main(int argc, char *argv[])
{

    ::ATORS::tendon_driven_actuation::ActuatedTendons actuated_tendons;

    const double d = 0.05;
    ::ATORS::tendon_driven_actuation::Tendon tendon1(Eigen::Vector3d(0, d, 0));
    ::ATORS::tendon_driven_actuation::Tendon tendon2(Eigen::Vector3d(0, -d, 0));


    const double C = 4;
    const double t_end = 10;

    const double omega = 2*M_PI*C/t_end;
    const double K1 = 30;
    const double K2 = 15;

    ::ATORS::simple_actuator::SimpleActuator motor_tendon_1([&](const double &t)->Eigen::VectorXd{
        const double theta = omega*t + M_PI;

        ::LieAlgebra::Vector1d tau_1;
        tau_1 << K1 * (1 + cos(theta));

        return tau_1;
    });

    ::ATORS::simple_actuator::SimpleActuator motor_tendon_2([&](const double &t)->Eigen::VectorXd{
        const double theta = omega*t + M_PI;

        ::LieAlgebra::Vector1d tau_2;
        tau_2 << K2 * (1 + cos(theta));
        return tau_2;
    });


    actuated_tendons = {
        { motor_tendon_1, tendon1 },
        { motor_tendon_2, tendon2 }
    };




    ::CROSP::InternallyActuatedCosseratRod rod(actuated_tendons);


    const unsigned int ne = rod.getCoordinatesDimension();
    Eigen::VectorXd q      = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd dot_q  = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

    rod.updateParameterisation(q, dot_q, ddot_q);

    rod.forwardKinematics();

    const auto Lambda_X1 = ::LieAlgebra::Vector6d::Zero();
    rod.backwardDynamics(Lambda_X1);

    rod.updateInternalActuation(1);

    const auto Q_ad = rod.getQad();

    std::cout << "Q_ad : \n" << Q_ad << "\n\n";

    return 0;
}
