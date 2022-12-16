#include "CROSP/CROSP/cosserat_rod.hpp"

#include <QApplication>
#include <QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>


#include "StaticInitialiser/static_initialiser.hpp"


#include "NEMI/newmark_integrator.hpp"

constexpr std::array<bool, 6> admitted_deformations = {
    false,
    true,
    false,

    false,
    false,
    false
};

constexpr unsigned int ne = 8;

constexpr unsigned int Nc = 21;



QApplication* app;
::Chebyshev::ChebyshevReconstructor rod_shape_reconstructor(Nc);

const double t_end = 0.2;
const double dt = 0.005;
::NewmarkIntegrator newmark_integrator(dt, t_end);





void printVector(const std::vector<double> &t_v,
                 const std::string &t_name)
{
    std::cout << t_name << " : ";
    for(const auto &e : t_v)
        std::cout << e << ", ";


    std::cout << "\n\n";
}


void findStaticSolution(::CROSP::CosseratRod &t_rod,
                        Eigen::VectorXd &t_q,
                        const ::LieAlgebra::Vector6d &t_Lambda_X1)
{

    t_rod.updateParameterisation(t_q, 0*t_q, 0*t_q);

    t_rod.forwardKinematics();

    t_rod.backwardDynamics(t_Lambda_X1);


    auto residual = t_rod.getStaticEquilibrium(t_q);


    Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(t_q.size());

    Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(t_q.size(), t_q.size());

    while(residual.norm() > 1e-6){

        for(unsigned int i=0; i<t_q.size(); i++){
            Delta_q.setZero();
            Delta_q[i] = 1;

            t_rod.updateParameterisationVariation(Delta_q, 0*Delta_q, 0*Delta_q);

            t_rod.forwardTangentKinematics();

            t_rod.backwardTangentDynamics(::LieAlgebra::Vector6d::Zero());


            Jacobian.col(i) = t_rod.getTangentStaticEquilibrium(Delta_q, 0*Delta_q);

        }

        t_q -= Jacobian.inverse() * residual;


        t_rod.updateParameterisation(t_q, 0*t_q, 0*t_q);

        t_rod.forwardKinematics();

        t_rod.backwardDynamics(t_Lambda_X1);


        residual = t_rod.getStaticEquilibrium(t_q);

    }

}




void plot(const ::CROSP::CosseratRod &t_rod)
{
    const Eigen::MatrixXd stack_of_points_rod =
            rod_shape_reconstructor.reconstructRodShape(t_rod.getRodPositionsAtChebyshevPoints());


    QtCharts::QChart* chart = new QtCharts::QChart();

    QtCharts::QLineSeries* series_rod = new QtCharts::QLineSeries();
    for(unsigned int j=0; j<stack_of_points_rod.rows(); j++)
        series_rod->append(stack_of_points_rod(j, 0), stack_of_points_rod(j, 2));
    chart->addSeries(series_rod);



    chart->createDefaultAxes();
    chart->axes(Qt::Horizontal).back()->setTitleText( "y [m]" );
    chart->axes(Qt::Horizontal).back()->setRange(-0.1, 1.1);
    chart->axes(Qt::Vertical).back()->setTitleText( "z [m]" );
    chart->axes(Qt::Vertical).back()->setRange(-0.4, 0.8);
    chart->legend();
    chart->setTitle( "Robot configuration" );

    QtCharts::QChartView* chart_view = new QtCharts::QChartView(chart);
    chart_view->resize(800,800);
    chart_view->show();

    app->exec();
}





int main(int argc, char *argv[])
{
    app = new QApplication(argc, argv);


    ::CROSP::polynomial_representation::PolynomialRepresentation poly(admitted_deformations, ne);



    const double mu = 0.0;
    const Eigen::Vector3d gravity = Eigen::Vector3d::Zero();
    ::CROSP::rod_properties::RodProperties rod_properties(mu, gravity, poly);




    ::CROSP::CosseratRod rod(poly, rod_properties, Nc);



    ::LieAlgebra::Vector6d F1;
    F1 <<   0,
            0,
            0,

            0,
            0,
            0.5;

    Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);


    findStaticSolution(rod, q, F1);
//    rod.updateParameterisation(q, 0*q, 0*q);
//    rod.forwardKinematics();
//    rod.backwardDynamics(F1);


    std::cout << "q : \n" << q << "\n\n";


    //plot(rod);


    const auto Lambda_stack = rod.m_idm_integrators->m_internal_couples->getStackAsMatrix();
    const auto Cy_stack = Lambda_stack(1, Eigen::all);



    std::cout << "Cy_stack : \n" << Cy_stack << "\n\n";

    const auto EI = rod.m_rod_properties->m_material_properties.m_E
                     * rod.m_rod_properties->m_rod_dimensions.m_Jy;


    Eigen::MatrixXd K_stack(3, Nc);
    for(unsigned int i=0; i<Nc; i++)
        K_stack.col(i) = rod.m_strain_parameterisation->m_K_stack->at((Nc-1)-i);


    std::cout << "K_stack : \n" << K_stack << "\n\n";

    const auto deformation = EI * K_stack.row(1);

    std::cout << "deformation : \n" << deformation << "\n\n";


    const auto error = Cy_stack - deformation;

    std::cout << "error : \n" << error << "\n\n";
















    //  Starting dynamic simulation

    std::vector<double> error_norm_stack;

    Eigen::VectorXd Delta_q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd Delta_dot_q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd Delta_ddot_q = Eigen::VectorXd::Zero(ne);

    Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(ne, ne);

    Eigen::VectorXd q_update = Eigen::VectorXd::Zero(ne);

    for(const auto simulation_time : newmark_integrator.getSimulationTimeSequence()){

        std::cout << "simulation time : " << simulation_time << std::endl;

        newmark_integrator.computeBallisticPrection(q, dot_q, ddot_q);

        rod.updateParameterisation(q, dot_q, ddot_q);

        rod.forwardKinematics();

        rod.backwardDynamics(::LieAlgebra::Vector6d::Zero());


        auto residual = rod.getStaticEquilibrium(q, dot_q);




        while(residual.norm() > 1e-6){

            for(unsigned int i=0; i<ne; i++){
                Delta_q.setZero();
                Delta_q[i] = 1;

                newmark_integrator.computeTangentCorrection(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.updateParameterisationVariation(Delta_q, Delta_dot_q, Delta_ddot_q);

                rod.forwardTangentKinematics();

                rod.backwardTangentDynamics(::LieAlgebra::Vector6d::Zero());


                Jacobian.col(i) = rod.getTangentStaticEquilibrium(Delta_q, Delta_dot_q);

            }

            q_update = -Jacobian.inverse() * residual;

            newmark_integrator.computeCorrection(q_update, q, dot_q, ddot_q);

            rod.updateParameterisation(q, dot_q, ddot_q);
            rod.forwardKinematics();
            rod.backwardDynamics(::LieAlgebra::Vector6d::Zero());


            residual = rod.getStaticEquilibrium(q, dot_q);

        }


        //  Solution found !!!


        const auto Lambda_stack = rod.m_idm_integrators->m_internal_couples->getStackAsMatrix();
        const auto Cy_stack = Lambda_stack(1, Eigen::all);


        Eigen::MatrixXd K_stack(3, Nc);
        for(unsigned int i=0; i<Nc; i++)
            K_stack.col(i) = rod.m_strain_parameterisation->m_K_stack->at((Nc-1)-i);


        const auto M_ang = rod.m_rod_properties->getMAngular();
        const auto dot_Omega_stack = rod.m_idm_integrators->m_angular_acceleration->getStackAsMatrix();


        const auto M_lin = rod.m_rod_properties->getMAngular();
        const auto V_stack = rod.m_idm_integrators->m_linear_velocity->getStackAsMatrix();


        Eigen::MatrixXd linear_acceleration_term(3, Nc);
        for(unsigned int j=0; j<Nc; j++){
            const auto V = V_stack.col(j);

            linear_acceleration_term.col(j) = ::LieAlgebra::skew( V ) * M_lin * V;
        }


        const Eigen::VectorXd inertial_term = (M_ang * dot_Omega_stack + linear_acceleration_term)(1, Eigen::all);

        const auto deformation = EI * K_stack.row(1);


        const auto error = Cy_stack - deformation - inertial_term;

        error_norm_stack.push_back( error.norm() );


    }



    printVector(error_norm_stack, "error_norm_stack");


    return 0;
}
