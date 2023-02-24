#include "CROSP/CROSP/cosserat_rod.hpp"














int main(int argc, char *argv[])
{


    const auto strain_par = std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>();


    ::CROSP::rod_properties::MaterialProperties mat_prop;
    ::CROSP::rod_properties::RodDimensions rod_dim;
    ::CROSP::rod_properties::RodProperties rod_properties(strain_par, rod_dim, mat_prop);


    ::CROSP::CosseratRod rod(strain_par, rod_properties);


    ::CROSP::CosseratRod rod2;

    const unsigned int ne = rod.getCoordinatesDimension();
    Eigen::VectorXd q      = Eigen::VectorXd::Random(ne);
    Eigen::VectorXd dot_q  = Eigen::VectorXd::Random(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Random(ne);

    rod.updateParameterisation(q, dot_q, ddot_q);

    rod.forwardKinematics();


    return 0;
}
