#include "strain_parameterisation.hpp"
#include "idm_integrators.hpp"
#include "cosserat_rod.hpp"
#include "tidm_integrators.hpp"



using namespace CROSP;




int main(int argc, char *argv[])
{

    std::shared_ptr<StrainParameterisation> strain_parameterisation = std::make_shared<StrainParameterisation>();

    Eigen::Matrix3d M_angular, M_linear;
    M_angular.setIdentity();
    M_linear.setIdentity();

    std::shared_ptr<IDMIntegrators> idm_integrators = std::make_shared<IDMIntegrators>(strain_parameterisation, M_angular, M_linear);

    return 0;
}
