
#include <iostream>

#include "cosserat_rod.hpp"
#include "strain_parameterisation.hpp"




using namespace CROSP;








int main(int argc, char *argv[])
{
//    const unsigned int number_of_Chebyshev_points = 11;

//    const auto Phi_stack = generatePhiStack(3, 3, ::Chebyshev::ComputeChebyshevPoints(number_of_Chebyshev_points));
//    StrainParameterisation strain_parameterisation(Phi_stack);

//    CosseratRod cosserat_rod(number_of_Chebyshev_points);

//    Eigen::VectorXd qe = Eigen::VectorXd::Zero(3*3);
//    Eigen::VectorXd dot_qe = Eigen::VectorXd::Zero(3*3);
//    Eigen::VectorXd ddot_qe = Eigen::VectorXd::Zero(3*3);

//    qe.setRandom();
//    dot_qe.setRandom();
//    ddot_qe.setRandom();

//    strain_parameterisation.update(qe, dot_qe, ddot_qe);

//    cosserat_rod.updateParameterisation(qe, dot_qe, ddot_qe);

//    std::cout << "Stacks : \n";
//    for(unsigned int i=0; i<number_of_Chebyshev_points; i++){
//        std::cout << "At point " << i << " the stakcs : \n";
//        std::cout << "  - K : \n";
//        std::cout << "          " << cosserat_rod.m_K_stack->at(i).transpose() << "\n";
//        std::cout << "          " << strain_parameterisation.m_K_stack->at(i).transpose() << "\n";
//        std::cout << "  - dot K : \n";
//        std::cout << "              " << cosserat_rod.m_dot_K_stack->at(i).transpose() << "\n";
//        std::cout << "              " << strain_parameterisation.m_dot_K_stack->at(i).transpose() << "\n";
//        std::cout << "  - ddot K : \n";
//        std::cout << "               " << cosserat_rod.m_ddot_K_stack->at(i).transpose() << "\n";
//        std::cout << "               " << strain_parameterisation.m_ddot_K_stack->at(i).transpose() << "\n";


//        std::cout << "\n\n";
//    }

    return 0;
}
