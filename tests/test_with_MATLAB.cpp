#include "CROSP/CROSP/cosserat_rod.hpp"

#include "utilities/Eigen/eigen_io.hpp"







int main(int argc, char *argv[])
{
    const std::string path = "../../../MATLAB/test_different_modes/";

    std::array<bool, 6> admitted_deformations = {
            true,
            true,
            true,
            false,
            false,
            false
        };

    Eigen::VectorXd deformations_stack = Eigen::VectorXd::Zero(6);
    for(unsigned int i=0; i<admitted_deformations.size(); i++)
        deformations_stack[i] = static_cast<unsigned int>( admitted_deformations[i] );
    writeToFile("deformations_stack", deformations_stack, path);

    std::vector<unsigned int> number_of_modes_stack {
        2, 3, 5
    };

    Eigen::VectorXd ne_stack = Eigen::VectorXd::Zero(6);
    unsigned int j =0;
    for(unsigned int i=0; i<admitted_deformations.size(); i++)
        if(admitted_deformations[i])
            ne_stack[i] = number_of_modes_stack[j++];
    writeToFile("ne_stack", ne_stack, path);






    const unsigned int number_of_Chebyshev_points = 31;


    CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, number_of_modes_stack);

    ::CROSP::rod_properties::RodDimensions rod_dimension;
    rod_dimension.m_L = 2.41;
    ::CROSP::CosseratRod<::CROSP::numerical_integrators::explicit_methods::ExplicitIntegrator<>> rod(polynomial_representation, number_of_Chebyshev_points, rod_dimension);



    writeToFile("Kee", rod.m_Kee, path);
    writeToFile("Dee", rod.m_Dee, path);


    writeToFile("B", polynomial_representation.m_B, path);
    writeToFile("B_bar", polynomial_representation.m_Bbar, path);


    const auto ne = polynomial_representation.getCoordinatesDimension();
    Eigen::VectorXd q      = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd dot_q  = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);


//    for(unsigned int i=0; i<ne; i++){
//        q(i)      = 1.0*(ne - i);
//        dot_q(i)  = 0.5*(ne - i);
//        ddot_q(i) = 0.2*(ne - i);
//    }

    q.setRandom();
    dot_q.setRandom();
    ddot_q.setRandom();

//   q(1)      = 1;
//   dot_q(1)  = 400;
//   ddot_q(1) = 160000;


//    dot_q(0) = 1;



    std::cout << "q : \n" << q << "\n\n" "dot q : \n" << dot_q << "\n\n" "ddot q : \n" << ddot_q << "\n\n";




    writeToFile("q", q, path);
    writeToFile("dot_q", dot_q, path);
    writeToFile("ddot_q", ddot_q, path);


    rod.m_cosserat_rod_integrators->updateParameterisation(q, dot_q, ddot_q);


    rod.m_cosserat_rod_integrators->forwardKinematics();



    auto zeros_6x1 = ::LieAlgebra::Vector6d::Zero();
    rod.m_cosserat_rod_integrators->backwardDynamics(zeros_6x1);

    const auto IDM_ODE_states_stacks = rod.m_cosserat_rod_integrators->getIDMStatesObservations();




    writeToFile("Q_stack", IDM_ODE_states_stacks.orientation_stack, path);
    writeToFile("r_stack", IDM_ODE_states_stacks.r_stack, path);
    writeToFile("Omega_stack", IDM_ODE_states_stacks.Omega_stack, path);
    writeToFile("V_stack", IDM_ODE_states_stacks.V_stack, path);
    writeToFile("dot_Omega_stack", IDM_ODE_states_stacks.dot_Omega_stack, path);
    writeToFile("dot_V_stack", IDM_ODE_states_stacks.dot_V_stack, path);
    writeToFile("N_stack", IDM_ODE_states_stacks.N_stack, path);
    writeToFile("C_stack", IDM_ODE_states_stacks.C_stack, path);
    writeToFile("Qa_stack", IDM_ODE_states_stacks.Qa_stack, path);
    writeToFile("Qad_stack", IDM_ODE_states_stacks.Qad_stack, path);



    Eigen::VectorXd Delta_q      = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd Delta_dot_q  = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


    const unsigned int a = 400;
    const unsigned int b = 160000;
//    const unsigned int a = 0;
//    const unsigned int b = 160000;

    for(unsigned int i=0; i<ne; i++){

        std::cout.flush();
        Delta_q.setZero();
        Delta_q[i] = 1;

        Delta_dot_q  = a * Delta_q;
        Delta_ddot_q = b * Delta_q;


        rod.m_cosserat_rod_integrators->updateDeltaParameterisation(Delta_q, Delta_dot_q, Delta_ddot_q);

        rod.m_cosserat_rod_integrators->forwardTangentKinematics();

        rod.m_cosserat_rod_integrators->backwardTangentDynamics(zeros_6x1);

        const auto TIDM_ODE_states_stacks = rod.m_cosserat_rod_integrators->getTIDMStatesObservations();


        writeToFile("Delta_rotation_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.orientation_stack, path);
        writeToFile("Delta_position_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.r_stack, path);
        writeToFile("Delta_Omega_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.Omega_stack, path);
        writeToFile("Delta_V_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.V_stack, path);
        writeToFile("Delta_dot_Omega_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.dot_Omega_stack, path);
        writeToFile("Delta_dot_V_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.dot_V_stack, path);
        writeToFile("Delta_N_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.N_stack, path);
        writeToFile("Delta_C_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.C_stack, path);
        writeToFile("Delta_Qa_stack" "_Delta" + std::to_string(i), TIDM_ODE_states_stacks.Qa_stack, path);

    }


    std::cout << "Ok saved all the data\n\n\n";
    std::cout.flush();

    return 0;
}
