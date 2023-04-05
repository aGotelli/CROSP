#include "CROSP/CROSP/cosserat_rod.hpp"





#include "CROSP/CROSP/internal_actuation.hpp"








int main(int argc, char *argv[])
{

    ::CROSP::CosseratRod rod;


    const unsigned int ne = rod.getCoordinatesDimension();
    Eigen::VectorXd q      = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd dot_q  = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);

    rod.updateParameterisation(q, dot_q, ddot_q);

    rod.forwardKinematics();

    const auto Lambda_X1 = ::LieAlgebra::Vector6d::Zero();
    rod.backwardDynamics(Lambda_X1);

    rod.internal_actuator.setActuation(1);

    const auto Q_ad = rod.internal_actuator.computedGeneraliseInternalActuation();

    std::cout << "Q_ad : \n" << Q_ad << "\n\n";

    return 0;
}
