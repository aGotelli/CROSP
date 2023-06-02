#include "CROSP/CROSP/cosserat_rod.hpp"


#include "utilities/Eigen/eigen_io.hpp"







int main(int argc, char *argv[])
{
    const std::string path = "../../../MATLAB/test_runge_kutta/";

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
        5, 2, 3
    };

    Eigen::VectorXd ne_stack = Eigen::VectorXd::Zero(6);
    unsigned int j =0;
    for(unsigned int i=0; i<admitted_deformations.size(); i++)
        if(admitted_deformations[i])
            ne_stack[i] = number_of_modes_stack[j++];
    writeToFile("ne_stack", ne_stack, path);






    const unsigned int number_of_Chebyshev_points = 31;
    auto poly_with_admitted_def_and_modes =
            std::make_shared<::CROSP::polynomial_representation::PolynomialRepresentation>(admitted_deformations, number_of_modes_stack);


    auto strain_parameterisation =
            std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(poly_with_admitted_def_and_modes, number_of_Chebyshev_points);

    ::CROSP::rod_properties::CircularCrossSection cs;
    const double length = 1.5;
    ::CROSP::rod_properties::RodDimensions rod_dimensions(&cs, length);

    //  Now use it in rod properties
    auto rod_properties =
            std::make_shared<::CROSP::rod_properties::RodProperties>(strain_parameterisation,
                                                                     rod_dimensions,
                                                                     ::CROSP::rod_properties::MaterialProperties());




    //  Variables related the perturbation of the strain parameterisation
    std::shared_ptr<::CROSP::strain_parameterisation::StrainParameterisation> strain_parameterisation_Delta {
        std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(strain_parameterisation,
                                                                          ::LieAlgebra::Vector6d::Zero())
    };


    ::CROSP::numerical_integrators::runge_kutta::RungeKuttaIntegrator m_integrators(*strain_parameterisation,
                                                                               *strain_parameterisation_Delta,
                                                                               *strain_parameterisation->m_polynomial_representation,
                                                                               *rod_properties,
                                                                               number_of_Chebyshev_points);




    writeToFile("Kee", rod_properties->m_Kee, path);
    writeToFile("Dee", rod_properties->m_Dee, path);


    writeToFile("B", strain_parameterisation->m_polynomial_representation->m_B, path);
    writeToFile("B_bar", strain_parameterisation->m_polynomial_representation->m_Bbar, path);


    const auto ne = strain_parameterisation->m_polynomial_representation->getCoordinatesDimension();
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


    std::cout << "q : \n" << q << "\n\n" "dot q : \n" << dot_q << "\n\n" "ddot q : \n" << ddot_q << "\n\n";




    writeToFile("q", q, path);
    writeToFile("dot_q", dot_q, path);
    writeToFile("ddot_q", ddot_q, path);


    m_integrators.updateParameterisation(q, dot_q, ddot_q);


    m_integrators.forwardKinematics();



    auto forward_kinematics_stack = m_integrators.m_forward_kinematics_stack;


    auto zeros = ::LieAlgebra::Vector6d::Zero();
    m_integrators.backwardDynamics(zeros);

    auto Lambda_Qa_stack = m_integrators.m_backward_Lambda_Qa_stack;



    const auto rod_shape = m_integrators.getRodPositions();


    std::cout << "Rod shape : \n" << rod_shape << "\n\n";

    const auto Q_stack         = forward_kinematics_stack.block<4, number_of_Chebyshev_points>(0, 0);
    const auto r_stack         = forward_kinematics_stack.block<3, number_of_Chebyshev_points>(4, 0);
    const auto Omega_stack     = forward_kinematics_stack.block<3, number_of_Chebyshev_points>(7, 0);
    const auto V_stack         = forward_kinematics_stack.block<3, number_of_Chebyshev_points>(10, 0);
    const auto dot_Omega_stack = forward_kinematics_stack.block<3, number_of_Chebyshev_points>(13, 0);
    const auto dot_V_stack     = forward_kinematics_stack.block<3, number_of_Chebyshev_points>(16, 0);

    const auto N_stack  = Lambda_Qa_stack.block<3, number_of_Chebyshev_points>(0, 0);
    const auto C_stack  = Lambda_Qa_stack.block<3, number_of_Chebyshev_points>(3, 0);
    const auto Qa_stack = Lambda_Qa_stack.block(6, 0, ne, 1);





    writeToFile("Q_stack", Q_stack, path);
    writeToFile("r_stack", r_stack, path);
    writeToFile("Omega_stack", Omega_stack, path);
    writeToFile("V_stack", V_stack, path);
    writeToFile("dot_Omega_stack", dot_Omega_stack, path);
    writeToFile("dot_V_stack", dot_V_stack, path);
    writeToFile("N_stack", N_stack, path);
    writeToFile("C_stack", C_stack, path);
    writeToFile("Qa_stack", Qa_stack, path);



    Eigen::VectorXd Delta_q      = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd Delta_dot_q  = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);


    const unsigned int a = 400;
    const unsigned int b = 160000;

    for(unsigned int i=0; i<ne; i++){

        Delta_q.setZero();
        Delta_q[i] = 1;

        Delta_dot_q  = a * Delta_q;
        Delta_ddot_q = b * Delta_q;


        m_integrators.updateDeltaParameterisation(Delta_q, Delta_dot_q, Delta_ddot_q);

        m_integrators.forwardTangentKinematics();

        m_integrators.backwardTangentDynamics(zeros);



//        const auto Delta_rotation_stack = m_integrators.m_tidm_integrators->m_Delta_rotation->getStackAsMatrix();
//        const auto Delta_position_stack = m_integrators.m_tidm_integrators->m_Delta_position->getStackAsMatrix();
//        const auto Delta_Omega_stack = m_integrators.m_tidm_integrators->m_Delta_angular_velocity->getStackAsMatrix();
//        const auto Delta_V_stack     = m_integrators.m_tidm_integrators->m_Delta_linear_velocity->getStackAsMatrix();
//        const auto Delta_dot_Omega_stack = m_integrators.m_tidm_integrators->m_Delta_angular_acceleration->getStackAsMatrix();
//        const auto Delta_dot_V_stack     = m_integrators.m_tidm_integrators->m_Delta_linear_acceleration->getStackAsMatrix();

//        const auto Delta_N_stack = m_integrators->m_tidm_integrators->m_Delta_internal_forces->getStackAsMatrix();
//        const auto Delta_C_stack = m_integrators->m_tidm_integrators->m_Delta_internal_couples->getStackAsMatrix();
//        const auto Delta_Qa_stack = m_integrators->m_tidm_integrators->m_Delta_generalised_forces->getStackAsMatrix();


//        writeToFile("Delta_rotation_stack" "_Delta" + std::to_string(i), Delta_rotation_stack, path);
//        writeToFile("Delta_position_stack" "_Delta" + std::to_string(i), Delta_position_stack, path);
//        writeToFile("Delta_Omega_stack" "_Delta" + std::to_string(i), Delta_Omega_stack, path);
//        writeToFile("Delta_V_stack" "_Delta" + std::to_string(i), Delta_V_stack, path);
//        writeToFile("Delta_dot_Omega_stack" "_Delta" + std::to_string(i), Delta_dot_Omega_stack, path);
//        writeToFile("Delta_dot_V_stack" "_Delta" + std::to_string(i), Delta_dot_V_stack, path);
//        writeToFile("Delta_N_stack" "_Delta" + std::to_string(i), Delta_N_stack, path);
//        writeToFile("Delta_C_stack" "_Delta" + std::to_string(i), Delta_C_stack, path);
//        writeToFile("Delta_Qa_stack" "_Delta" + std::to_string(i), Delta_Qa_stack, path);

    }



    std::cout << "Ok saved all the data\n\n\n";
    std::cout.flush();

    return 0;
}
