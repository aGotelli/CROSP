
#include "cosserat_rod.hpp"




int main(int argc, char *argv[])
{

    CROSP::CosseratRod cosserat_rod;

    constexpr unsigned int ne = 3;
    constexpr unsigned int na = 3;
    Eigen::Matrix<double, ne*na, 1> qe;
    //  Here we give some value for the strain
    qe <<   0,
            0,
            0,
            1.2877691307032,
           -1.63807499160786,
            0.437406679142598,
            0,
            0,
            0;

    Eigen::Matrix<double, ne*na, 1> dot_qe, ddot_qe;
    dot_qe.setZero();
    ddot_qe.setZero();

    cosserat_rod.updateParameterisation(qe, dot_qe, ddot_qe);

    cosserat_rod.forwardKinematics();



    return 0;
}
