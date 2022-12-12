#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include <benchmark/benchmark.h>



constexpr unsigned int repetitions = 20;

constexpr std::array<bool, 6> admitted_deformations = {
    false,
    true,
    false,

    false,
    false,
    false
};

constexpr unsigned int ne = 4;



void solutionSpectral()
{
    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation);


    Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
    q[0] = 1;
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);


    rod.updateParameterisation(q, dot_q, ddot_q);

    rod.m_idm_integrators->m_quaternion->solveSystem();
    const auto Q_X1 = rod.m_idm_integrators->m_quaternion->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);
    const auto R_spectral = Eigen::Quaterniond(Q_X1(0),
                                               Q_X1(1),
                                               Q_X1(2),
                                               Q_X1(3)).toRotationMatrix();

    std::cout << "R(X=1) spectral \n" << R_spectral << std::endl;

}


void compareSolutions()
{

    const unsigned int Nc = 25;

    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation, Nc);


    Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
    q[0] = 1.0;
    q[1] = 0.8;
    q[2] = -0.3;
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);


    rod.updateParameterisation(q, dot_q, ddot_q);

    rod.m_idm_integrators->m_quaternion->solveSystem();
    const auto Q_X1 = rod.m_idm_integrators->m_quaternion->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);
    const auto R_spectral = Eigen::Quaterniond(Q_X1(0),
                                               Q_X1(1),
                                               Q_X1(2),
                                               Q_X1(3)).toRotationMatrix();




    const double number_of_Chebyshev_points = rod.m_strain_parameterisation->m_number_of_Chebyshev_points;
    const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(number_of_Chebyshev_points);


    double begin, end;
    double h, step;

    double x1, x2, x3;

    ::LieAlgebra::Vector6d Xi1;
    ::LieAlgebra::Vector6d Xi2;
    ::LieAlgebra::Vector6d Xi3;

    Eigen::Matrix3d A1;
    Eigen::Matrix3d A2;
    Eigen::Matrix3d A3;

    Eigen::Matrix3d alpha1;
    Eigen::Matrix3d alpha2;
    Eigen::Matrix3d alpha3;

    Eigen::Matrix3d C1;
    Eigen::Matrix3d C2;

    Eigen::Matrix3d Omega;

    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();


    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){


        begin = Chebyshev_points[ (number_of_Chebyshev_points-1) - i];
        end = Chebyshev_points[ (number_of_Chebyshev_points-1) - (i + 1)];

        h = end - begin;
        step = h*sqrt(15.0)/10.0;

        x1 = begin + (0.5 - step)*h;
        x2 = begin + 0.5*h;
        x3 = begin + (0.5 + step)*h;


        Xi1 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1) * q;
        Xi2 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2) * q;
        Xi3 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3) * q;


        A1 = ::LieAlgebra::skew( Xi1.block<3,1>(0, 0) );
        A2 = ::LieAlgebra::skew( Xi2.block<3,1>(0, 0) );
        A3 = ::LieAlgebra::skew( Xi3.block<3,1>(0, 0) );


        alpha1 = h*A2;
        alpha2 = (A3 - A1)*sqrt(15.0)*h/3.0;
        alpha3 = (A3 - 2*A2 + A1)*10.0*h/3.0;


        C1 = alpha1*alpha2 - alpha2*alpha1;
        C2 = -( alpha1*(2*alpha3 + C1) - (2*alpha3 + C1)*alpha1 )/60.0;

        Omega = alpha1
                + alpha3/12.0
                + ( (-20*alpha1 - alpha3 + C1)*(alpha2 + C2)
                    - (alpha2 + C2)*(-20*alpha1 - alpha3 + C1) )/240.0;


        R = R * ::LieAlgebra::expRodigues(Omega);

    }




    const auto error = ::LieAlgebra::logSO3(R, R_spectral);


    //  Compare the solutions
    std::cout << "With q : \n" << q << "\n" <<
                 "R(X=1) OSNI   : \n" << R_spectral << "\n" <<
                 "R(X=1) Magnus : \n" << R << "\n" <<
                 "\n" <<
                 "error : \n" << error << "\n" <<
                 std::endl;

}





int main(int argc, char *argv[])
{

    compareSolutions();



    Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
    q[0] = 1;
    Eigen::VectorXd dot_q = Eigen::VectorXd::Zero(ne);
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Zero(ne);




    ::benchmark::RegisterBenchmark("Quaternion spectral", [&](::benchmark::State &t_state) {


        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        rod.updateParameterisation(q, dot_q, ddot_q);

        while(t_state.KeepRunning()){
            rod.m_idm_integrators->m_quaternion->solveSystem();

            if(std::isnan(rod.m_idm_integrators->m_quaternion->getStateAtPoint(0).x()) ||
                std::isnan(rod.m_idm_integrators->m_quaternion->getStateAtPoint(0).y()) ||
                std::isnan(rod.m_idm_integrators->m_quaternion->getStateAtPoint(0).z()) ||
                std::isnan(rod.m_idm_integrators->m_quaternion->getStateAtPoint(0).w()) )
                t_state.SkipWithError("Result is nan!");
        }


    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);




    ::benchmark::RegisterBenchmark("Quaternion spectral", [&](::benchmark::State &t_state) {

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);


        const double number_of_Chebyshev_points = rod.m_strain_parameterisation->m_number_of_Chebyshev_points;
        const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(number_of_Chebyshev_points);


        double begin, end;
        double h, step;

        double x1, x2, x3;

        ::LieAlgebra::Vector6d Xi1;
        ::LieAlgebra::Vector6d Xi2;
        ::LieAlgebra::Vector6d Xi3;

        Eigen::Matrix3d A1;
        Eigen::Matrix3d A2;
        Eigen::Matrix3d A3;

        Eigen::Matrix3d alpha1;
        Eigen::Matrix3d alpha2;
        Eigen::Matrix3d alpha3;

        Eigen::Matrix3d C1;
        Eigen::Matrix3d C2;

        Eigen::Matrix3d Omega;

        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();


        while(t_state.KeepRunning()){
            for(unsigned int i=0; i<number_of_Chebyshev_points -1; i++){


                begin = Chebyshev_points[ (number_of_Chebyshev_points-1) - i];
                end = Chebyshev_points[ (number_of_Chebyshev_points-1) - (i + 1)];

                h = end - begin;
                step = h*sqrt(15.0)/10.0;

                x1 = begin + (0.5 - step)*h;
                x2 = begin + 0.5*h;
                x3 = begin + (0.5 + step)*h;


                Xi1 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1) * q;
                Xi2 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2) * q;
                Xi3 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3) * q;


                A1 = ::LieAlgebra::skew( Xi1.block<3,1>(0, 0) );
                A2 = ::LieAlgebra::skew( Xi2.block<3,1>(0, 0) );
                A3 = ::LieAlgebra::skew( Xi3.block<3,1>(0, 0) );


                alpha1 = h*A2;
                alpha2 = (A3 - A1)*sqrt(15.0)*h/3.0;
                alpha3 = (A3 - 2*A2 + A1)*10.0*h/3.0;


                C1 = alpha1*alpha2 - alpha2*alpha1;
                C2 = ( alpha1*(2*alpha3 + C1) - (2*alpha3 + C1)*alpha1 );

                Omega = alpha1
                        + alpha3/12
                        + ( (-20*alpha1 - alpha3 + C1)*(alpha2 + C2)
                            - (alpha2 + C2)*(-20*alpha1 - alpha3 + C1) )/240;


                R = R * ::LieAlgebra::expRodigues(Omega);

            }

            if(R.hasNaN())
                t_state.SkipWithError("Result is nan!");
        }


//        const auto error = ::LieAlgebra::logSO3(R, R_spectral);

//        t_state.counters["norm(logSO(3))"] = error.norm();

        //std::cout << "R(X=1) magnus \n" << R << std::endl;
    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);




//    ::benchmark::Initialize(&argc, argv);


//    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
