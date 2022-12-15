#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include <benchmark/benchmark.h>


#include <boost/numeric/odeint.hpp>
#include <boost/range/algorithm.hpp>



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



Eigen::MatrixXd solutionSpectral(::CROSP::CosseratRod t_rod,
                                 const Eigen::VectorXd &t_q)
{


    t_rod.updateParameterisation(t_q, 0*t_q, 0*t_q);

    t_rod.m_idm_integrators->m_quaternion->solveSystem();
    const auto Q_X1 = t_rod.m_idm_integrators->m_quaternion->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);
    const auto R = Eigen::Quaterniond(Q_X1(0),
                                        Q_X1(1),
                                        Q_X1(2),
                                        Q_X1(3)).toRotationMatrix();

    return R;

}


Eigen::Matrix4d getA(const Eigen::Vector3d t_k)
{
    Eigen::Matrix4d A;
    A   <<     0  , -t_k(0),  -t_k(1),  -t_k(2),
            t_k(0),     0  ,   t_k(2),  -t_k(1),
            t_k(1), -t_k(2),      0  ,   t_k(0),
            t_k(2),  t_k(1),  -t_k(0),      0  ;

    return A;
}


Eigen::Matrix3d solutionODE(const ::CROSP::CosseratRod &rod,
                            const Eigen::VectorXd &q)
{

    typedef Eigen::Matrix<double, 4, 1 > state_type;

    typedef boost::numeric::odeint::runge_kutta_dopri5< state_type, double,
                                                        state_type, double,
                                                        boost::numeric::odeint::vector_space_algebra> stepper;


    state_type y;
    y << 1, 0, 0, 0;



    ::LieAlgebra::Vector6d Xi;
    Eigen::Vector3d k;


    const double L = 1.0;
    const double dX = 0.005;

    boost::numeric::odeint::integrate_adaptive(stepper(),
                                            [&](const state_type& t_y, state_type& t_dyds, const double t_X){

        Xi = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(t_X) * q;

        k = Xi.block<3,1>(0,0);

        t_dyds = 0.5*getA(k)*t_y;

    }, y, 0.0, L, dX);


    const Eigen::Quaterniond Q(y[0], y[1],y[2], y[3]);


    return Q.toRotationMatrix();

}



Eigen::MatrixXd solutionMagnus(const ::CROSP::CosseratRod &rod,
                               const Eigen::VectorXd &q)
{


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
        step = sqrt(15.0)/10.0;

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

    return R;

}



Eigen::MatrixXd solutionCayley(const ::CROSP::CosseratRod &rod,
                               const Eigen::VectorXd &q)
{


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
    Eigen::Matrix3d Omega_square;

    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d C;

    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();


    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){


        begin = Chebyshev_points[ (number_of_Chebyshev_points-1) - i];
        end = Chebyshev_points[ (number_of_Chebyshev_points-1) - (i + 1)];

        h = end - begin;
        step = sqrt(15.0)/10.0;

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

        Omega_square = Omega*Omega;
        C = Omega*(I - Omega_square*(I - Omega_square/10.0)/12.0);

        R = R * (I - C/2.0).inverse() * (I + C/2.0);

    }

    return R;

}


void compareSolutions()
{


    const unsigned int Nc = 31;

    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation, Nc);





    Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
    q.setRandom();


    const auto R_ODE    = solutionODE(rod, q);
    const auto R_OSNI   = solutionSpectral(rod, q);
    const auto R_Magnus = solutionMagnus(rod, q);
    const auto R_Cayley = solutionCayley(rod, q);


    std::cout << "With q = [" << q.transpose() << "]\n\n";


    std::cout << "Solutions :\n" <<
                 "  OSNI\n" <<
                 "      R : \n" << R_OSNI << "\n\n" <<
                 "  ODE\n" <<
                 "      R : \n" << R_ODE << "\n\n" <<
                 "  Magnus\n" <<
                 "      R : \n" << R_Magnus << "\n\n" <<
                 "  Cayley\n" <<
                 "      R : \n" << R_Cayley << "\n\n";

    std::cout << "\n\n\n";


    Eigen::Matrix3d RTR_Magnus = R_Magnus.transpose()*R_Magnus;

    std::cout << "R^T R :\n"<<
                 "  OSNI\n" <<
                 R_OSNI.transpose()*R_OSNI << "\n\n" <<
                 "  ODE\n" <<
                 R_ODE.transpose()*R_ODE << "\n\n" <<
                 "  Magnus\n" <<
                 R_Magnus.transpose()*R_Magnus << "\n\n" <<
                 "  Cayley\n" <<
                 R_Cayley.transpose()*R_Cayley << "\n\n";

    const auto error_ODE_OSNI = ::LieAlgebra::logSO3(R_ODE, R_OSNI);
    const auto error_ODE_Magnus = ::LieAlgebra::logSO3(R_ODE, R_Magnus);
    const auto error_ODE_Cayley = ::LieAlgebra::logSO3(R_ODE, R_Cayley);

    std::cout << "\n\n\n";

    std::cout << "Errors :\n" <<
                 "  ODE-OSNI\n" <<
                 error_ODE_OSNI << "\n\n" <<
                 "  ODE-Magnus\n" <<
                 error_ODE_Magnus << "\n\n" <<
                 "  ODE-Cayley\n" <<
                 error_ODE_Cayley << "\n\n";

}





int main(int argc, char *argv[])
{

    compareSolutions();


    ::benchmark::RegisterBenchmark("Quaternion ODE", [&](::benchmark::State &t_state) {

        const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        rod.updateParameterisation(q, 0*q, 0*q);


        typedef Eigen::Matrix<double, 4, 1 > state_type;

        typedef boost::numeric::odeint::runge_kutta_dopri5< state_type, double,
                                                            state_type, double,
                                                            boost::numeric::odeint::vector_space_algebra> stepper;


        state_type y;
        y << 1, 0, 0, 0;



        ::LieAlgebra::Vector6d Xi;
        Eigen::Vector3d k;


        const double L = 1.0;
        const double dX = 0.005;

        while(t_state.KeepRunning()){
            boost::numeric::odeint::integrate_adaptive(stepper(),
                                                    [&](const state_type& t_y, state_type& t_dyds, const double t_X){

                Xi = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(t_X) * q;

                k = Xi.block<3,1>(0,0);

                t_dyds = 0.5*getA(k)*t_y;

            }, y, 0.0, L, dX);
        }


    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);


    ::benchmark::RegisterBenchmark("Quaternion OSNI", [&](::benchmark::State &t_state) {

        const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        rod.updateParameterisation(q, 0*q, 0*q);

        while(t_state.KeepRunning()){
            rod.m_idm_integrators->m_quaternion->solveSystem();
        }


    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);




    ::benchmark::RegisterBenchmark("Quaternion Magnus", [&](::benchmark::State &t_state) {

        const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        rod.updateParameterisation(q, 0*q, 0*q);


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
        }

    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);




    ::benchmark::RegisterBenchmark("Quaternion Cayley", [&](::benchmark::State &t_state) {

        const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation);

        rod.updateParameterisation(q, 0*q, 0*q);


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
        Eigen::Matrix3d Omega_square;

        const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d C;

        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();


        while(t_state.KeepRunning()){
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


                Omega_square = Omega*Omega;
                C = Omega*(I - Omega_square*(I - Omega_square/10.0)/12.0);

                R = R * (I - C/2.0).inverse() * (I + C/2.0);

            }
        }

    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
