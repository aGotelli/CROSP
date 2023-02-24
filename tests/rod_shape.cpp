
#include "CROSP/CROSP/cosserat_rod.hpp"

#include <QApplication>
#include <QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>


#include "utilities/Eigen/eigen_io.hpp"


void plot(const ::CROSP::CosseratRod &t_rod);

const unsigned int number_of_Chebyshev_points = 31;

QApplication* app;
::Chebyshev::ChebyshevReconstructor rod_shape_reconstructor =
        ::Chebyshev::ChebyshevReconstructor(number_of_Chebyshev_points);

int main(int argc, char *argv[])
{

    app = new QApplication(argc, argv);

    const unsigned int ne = 4;


    ::CROSP::CosseratRod rod(number_of_Chebyshev_points);


    Eigen::VectorXd qe(ne);
    qe <<  2,
          0,
          1,
          0;

    rod.updateParameterisation(qe, 0*qe, 0*qe);
    rod.forwardKinematics();

    plot(rod);


    Eigen::MatrixXd rod_shape_reconstructed(rod_shape_reconstructor.getNumberOfReconstructionPoints(), 3);

    rod_shape_reconstructed = rod_shape_reconstructor.reconstructRodShape( rod.getRodPositionsAtChebyshevPoints() );




    writeToFile("rod_shape", rod_shape_reconstructed, "../../rod_shape/");



    return 0;
}


void plot(const ::CROSP::CosseratRod &t_rod)
{
    const Eigen::MatrixXd stack_of_points_rod =
            rod_shape_reconstructor.reconstructRodShape(t_rod.getRodPositionsAtChebyshevPoints());




    QtCharts::QChart* chart = new QtCharts::QChart();

    QtCharts::QLineSeries* series_rod = new QtCharts::QLineSeries();
    for(unsigned int j=0; j<stack_of_points_rod.rows(); j++)
        series_rod->append(stack_of_points_rod(j, 0), stack_of_points_rod(j, 1));
    chart->addSeries(series_rod);



    chart->createDefaultAxes();
    chart->axes(Qt::Horizontal).back()->setTitleText( "y [m]" );
    chart->axes(Qt::Horizontal).back()->setRange(-0.1, 1.1);
    chart->axes(Qt::Vertical).back()->setTitleText( "z [m]" );
    chart->axes(Qt::Vertical).back()->setRange(-0.1, 1.1);
    chart->legend();
    chart->setTitle( "Robot configuration" );

    QtCharts::QChartView* chart_view = new QtCharts::QChartView(chart);
    chart_view->resize(800,800);
    chart_view->show();


    app->exec();
}
