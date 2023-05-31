#include "CROSP/CROSP/cosserat_rod.hpp"
#include "CROSP/CROSP/internally_actuated_cosserat_rod.hpp"

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


    auto strain_par =
            std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(poly_with_admitted_def_and_modes, number_of_Chebyshev_points);

    ::CROSP::rod_properties::CircularCrossSection cs;
    const double length = 1.5;
    ::CROSP::rod_properties::RodDimensions rod_dimensions(&cs, length);

    //  Now use it in rod properties
    auto rod_properties =
            std::make_shared<::CROSP::rod_properties::RodProperties>(strain_par,
                                                                     rod_dimensions,
                                                                     ::CROSP::rod_properties::MaterialProperties());




    ::CROSP::CosseratRod rod(strain_par, rod_properties);
    rod.printProperties();


    writeToFile("Kee", rod.m_rod_properties->m_Kee, path);
    writeToFile("Dee", rod.m_rod_properties->m_Dee, path);


    writeToFile("B", rod.m_strain_parameterisation->m_polynomial_representation->m_B, path);
    writeToFile("B_bar", rod.m_strain_parameterisation->m_polynomial_representation->m_Bbar, path);


    const auto ne = rod.getCoordinatesDimension();
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


//    rod.updateParameterisation(q, dot_q, ddot_q);
    rod.m_idm_integrators->m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);
    rod.m_tidm_integrators->m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);


    rod.forwardKinematics();

    auto zeros = ::LieAlgebra::Vector6d::Zero();
    rod.backwardDynamics(zeros);

    const auto Q_stack = rod.m_idm_integrators->m_quaternion->getStackAsMatrix();
    const auto r_stack = rod.m_idm_integrators->m_position->getStackAsMatrix();
    const auto Omega_stack = rod.m_idm_integrators->m_angular_velocity->getStackAsMatrix();
    const auto V_stack     = rod.m_idm_integrators->m_linear_velocity->getStackAsMatrix();
    const auto dot_Omega_stack = rod.m_idm_integrators->m_angular_acceleration->getStackAsMatrix();
    const auto dot_V_stack     = rod.m_idm_integrators->m_linear_acceleration->getStackAsMatrix();

    const auto N_stack = rod.m_idm_integrators->m_internal_forces->getStackAsMatrix();
    const auto C_stack = rod.m_idm_integrators->m_internal_couples->getStackAsMatrix();
    const auto Qa_stack = rod.m_idm_integrators->m_generalised_forces->getStackAsMatrix();





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


//        rod.updateDeltaParameterisation(Delta_q, Delta_dot_q, Delta_ddot_q);
        rod.m_tidm_integrators->m_Delta_strain_parameterisation_stack->updateStrainParameterisation(Delta_q, Delta_dot_q, Delta_ddot_q);

        rod.forwardTangentKinematics();

        rod.backwardTangentDynamics(zeros);



        const auto Delta_rotation_stack = rod.m_tidm_integrators->m_Delta_rotation->getStackAsMatrix();
        const auto Delta_position_stack = rod.m_tidm_integrators->m_Delta_position->getStackAsMatrix();
        const auto Delta_Omega_stack = rod.m_tidm_integrators->m_Delta_angular_velocity->getStackAsMatrix();
        const auto Delta_V_stack     = rod.m_tidm_integrators->m_Delta_linear_velocity->getStackAsMatrix();
        const auto Delta_dot_Omega_stack = rod.m_tidm_integrators->m_Delta_angular_acceleration->getStackAsMatrix();
        const auto Delta_dot_V_stack     = rod.m_tidm_integrators->m_Delta_linear_acceleration->getStackAsMatrix();

        const auto Delta_N_stack = rod.m_tidm_integrators->m_Delta_internal_forces->getStackAsMatrix();
        const auto Delta_C_stack = rod.m_tidm_integrators->m_Delta_internal_couples->getStackAsMatrix();
        const auto Delta_Qa_stack = rod.m_tidm_integrators->m_Delta_generalised_forces->getStackAsMatrix();


        writeToFile("Delta_rotation_stack" "_Delta" + std::to_string(i), Delta_rotation_stack, path);
        writeToFile("Delta_position_stack" "_Delta" + std::to_string(i), Delta_position_stack, path);
        writeToFile("Delta_Omega_stack" "_Delta" + std::to_string(i), Delta_Omega_stack, path);
        writeToFile("Delta_V_stack" "_Delta" + std::to_string(i), Delta_V_stack, path);
        writeToFile("Delta_dot_Omega_stack" "_Delta" + std::to_string(i), Delta_dot_Omega_stack, path);
        writeToFile("Delta_dot_V_stack" "_Delta" + std::to_string(i), Delta_dot_V_stack, path);
        writeToFile("Delta_N_stack" "_Delta" + std::to_string(i), Delta_N_stack, path);
        writeToFile("Delta_C_stack" "_Delta" + std::to_string(i), Delta_C_stack, path);
        writeToFile("Delta_Qa_stack" "_Delta" + std::to_string(i), Delta_Qa_stack, path);

    }

    ::ATORS::tendon_driven_actuation::ActuatedTendons actuated_tendons;

    const double d_tendon = -0.011;
    ::ATORS::tendon_driven_actuation::Tendon tendon1(Eigen::Vector3d(0, d_tendon, 0));


    unsigned int weight = 2000;
    const double gravity_Toronto = 9.80436;

    ::ATORS::simple_actuator::SimpleActuator motor_tendon_1([&](const double &t)->Eigen::VectorXd{
        ::LieAlgebra::Vector1d tension;
        tension << weight*gravity_Toronto/1000;

        return tension;
    });




    actuated_tendons = {
        { motor_tendon_1, tendon1 }
    };


    //  Now the internal actuation
    ::ATORS::tendon_driven_actuation::TendonDrivenActuation tendond_driven_actuation(strain_par->m_K_stack,
                                                                                    strain_par->m_Gamma_stack,
                                                                                    strain_par->m_map_to_strain_stack,
                                                                                    actuated_tendons,
                                                                                    rod_dimensions.m_L);

    tendond_driven_actuation.updateActuation(0.0);
    const auto Q_ad_stack = tendond_driven_actuation.getStackAsMatrix();
    Eigen::VectorXd Q_ad = tendond_driven_actuation.getActuation();

    std::cout << "Q_ad:\n" << Q_ad << "\n\n";

    writeToFile("Q_ad_stack" + std::to_string(weight) + "g", Q_ad_stack, path);


    std::cout << "Ok saved all the data\n\n\n";
    std::cout.flush();

    return 0;
}
