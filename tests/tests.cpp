
//#include "cosserat_rod.hpp"

#include <memory>

#include "OSNI/include/ode_solver_interface.hpp"
#include "OSNI/include/ode_a.hpp"
#include "OSNI/include/ode_ab.hpp"
#include "OSNI/include/ode_b.hpp"

struct ODE : public OSNI::ODEA {


    ODE() : OSNI::ODEA(4) {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)override
    {
        return Eigen::MatrixXd::Zero(4, 4);

    }
};


int main(int argc, char *argv[])
{

    //QuaternionIntegrator q_int;


    std::shared_ptr<OSNI::ODESolverInterface> ode_ptr;// = std::make_shared<OSNI::ODEA>( 4 );

    ODE ode;



    return 0;
}
