#include "CROSP/CROSP/cosserat_rod.hpp"

#include <QApplication>
#include <QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>


#include "StaticInitialiser/static_initialiser.hpp"

constexpr std::array<bool, 6> admitted_deformations = {
    false,
    true,
    false,

    false,
    false,
    false
};

constexpr unsigned int ne = 10;

constexpr unsigned int Nc = 21;



QApplication* app;
::Chebyshev::ChebyshevReconstructor rod_shape_reconstructor(Nc);








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


    return 0;
}
