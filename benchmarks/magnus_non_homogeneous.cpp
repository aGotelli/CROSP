#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include <benchmark/benchmark.h>


#include <boost/numeric/odeint.hpp>
#include <boost/range/algorithm.hpp>


#include <eigen3/unsupported/Eigen/MatrixFunctions>

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

constexpr unsigned int Nc = 17;



Eigen::Vector3d solutionSpectral(::CROSP::CosseratRod t_rod,
                                 const Eigen::VectorXd &t_q,
                                 const Eigen::VectorXd &t_dot_q)
{


    t_rod.updateParameterisation(t_q, t_dot_q, 0*t_dot_q);

    t_rod.m_idm_integrators->m_angular_velocity->solveSystem();
    const auto Omega_X1 = t_rod.m_idm_integrators->m_angular_velocity->getStateAtPoint(::OSNI::INTEGRATION_DOMAIN::END);

    return Omega_X1;

}





Eigen::Vector3d solutionMagnus(const ::CROSP::CosseratRod &rod,
                               const Eigen::VectorXd &t_q,
                               const Eigen::VectorXd &t_dot_q)
{


    const double number_of_Chebyshev_points = rod.m_strain_parameterisation->m_number_of_Chebyshev_points;
    const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(number_of_Chebyshev_points);


    double begin, end;
    double h, step;

    double x1, x2, x3;

    ::LieAlgebra::Vector6d Xi1;
    ::LieAlgebra::Vector6d Xi2;
    ::LieAlgebra::Vector6d Xi3;


    ::LieAlgebra::Vector6d dot_Xi1;
    ::LieAlgebra::Vector6d dot_Xi2;
    ::LieAlgebra::Vector6d dot_Xi3;


    Eigen::Matrix3d M1;
    Eigen::Matrix3d M2;
    Eigen::Matrix3d M3;

    Eigen::Vector3d h1;
    Eigen::Vector3d h2;
    Eigen::Vector3d h3;


    Eigen::Matrix4d A1 = Eigen::Matrix4d::Zero();
    Eigen::Matrix4d A2 = Eigen::Matrix4d::Zero();
    Eigen::Matrix4d A3 = Eigen::Matrix4d::Zero();

    Eigen::Matrix4d alpha1;
    Eigen::Matrix4d alpha2;
    Eigen::Matrix4d alpha3;

    Eigen::Matrix4d C1;
    Eigen::Matrix4d C2;

    Eigen::Matrix4d Omega;

    Eigen::Vector4d X(0, 0, 0, 1);


    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){


        begin = Chebyshev_points[ (number_of_Chebyshev_points-1) - i];
        end = Chebyshev_points[ (number_of_Chebyshev_points-1) - (i + 1)];

        h = end - begin;
        step = sqrt(15.0)/10.0;

        x1 = begin + (0.5 - step)*h;
        x2 = begin + 0.5*h;
        x3 = begin + (0.5 + step)*h;


        Xi1 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1) * t_q;
        Xi2 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2) * t_q;
        Xi3 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3) * t_q;


        M1 = ::LieAlgebra::skew( Xi1.block<3,1>(0, 0) );
        M2 = ::LieAlgebra::skew( Xi2.block<3,1>(0, 0) );
        M3 = ::LieAlgebra::skew( Xi3.block<3,1>(0, 0) );




        dot_Xi1 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1) * t_dot_q;
        dot_Xi2 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2) * t_dot_q;
        dot_Xi3 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3) * t_dot_q;

        h1 = dot_Xi1.block<3, 1>(0, 0);
        h2 = dot_Xi2.block<3, 1>(0, 0);
        h3 = dot_Xi3.block<3, 1>(0, 0);




        A1.block<3, 4>(0, 0) << M1, h1;
        A2.block<3, 4>(0, 0) << M2, h2;
        A3.block<3, 4>(0, 0) << M3, h3;



        alpha1 = h*A2;
        alpha2 = (A3 - A1)*sqrt(15.0)*h/3.0;
        alpha3 = (A3 - 2*A2 + A1)*10.0*h/3.0;


        C1 = alpha1*alpha2 - alpha2*alpha1;
        C2 = -( alpha1*(2*alpha3 + C1) - (2*alpha3 + C1)*alpha1 )/60.0;

        Omega = alpha1
                + alpha3/12.0
                + ( (-20*alpha1 - alpha3 + C1)*(alpha2 + C2)
                    - (alpha2 + C2)*(-20*alpha1 - alpha3 + C1) )/240.0;


        X = Omega.exp() * X;

    }




    return X.block<3, 1>(0, 0);

}




void compareSolutions()
{


    const ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation, Nc);




    const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);
    const Eigen::VectorXd dot_q = Eigen::VectorXd::Random(ne);


    const auto Omega_X1_OSNI   = solutionSpectral(rod, q, dot_q);
    const auto Omega_X1_Magnus = solutionMagnus(rod, q, dot_q);

    std::cout << "With \n" <<
                 "      q = [" << q.transpose() << "]\n" <<
                 "  dot_q = [" << dot_q.transpose() << "]\n\n";


    std::cout << "Solutions :\n" <<
                 "  OSNI\n" <<
                 "      Omega : \n" << Omega_X1_OSNI << "\n\n" <<
                 "  Magnus\n" <<
                 "      Omega : \n" << Omega_X1_Magnus << "\n\n";

    std::cout << "\n\n\n";

    const auto error_OSNI_Magnus = Omega_X1_Magnus - Omega_X1_OSNI;

    std::cout << "Errors :\n" <<
                 "  OSNI-Magnus\n" <<
                 error_OSNI_Magnus << "\n\n";


}





int main(int argc, char *argv[])
{


    srand(time(nullptr));

    compareSolutions();




    ::benchmark::RegisterBenchmark("Omega OSNI", [](::benchmark::State &t_state) {

        const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);
        const Eigen::VectorXd dot_q = Eigen::VectorXd::Random(ne);

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation, Nc);

        rod.updateParameterisation(q, dot_q, 0*dot_q);


        while(t_state.KeepRunning()){
            rod.m_idm_integrators->m_angular_velocity->solveSystem();
        }


    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);




    ::benchmark::RegisterBenchmark("Omega Magnus", [&](::benchmark::State &t_state) {

        const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);
        const Eigen::VectorXd dot_q = Eigen::VectorXd::Random(ne);

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation, Nc);

        rod.updateParameterisation(q, dot_q, 0*dot_q);



        const double number_of_Chebyshev_points = rod.m_strain_parameterisation->m_number_of_Chebyshev_points;
        const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(number_of_Chebyshev_points);


        double begin, end;
        double h, step;

        double x1, x2, x3;

        ::LieAlgebra::Vector6d Xi1;
        ::LieAlgebra::Vector6d Xi2;
        ::LieAlgebra::Vector6d Xi3;


        ::LieAlgebra::Vector6d dot_Xi1;
        ::LieAlgebra::Vector6d dot_Xi2;
        ::LieAlgebra::Vector6d dot_Xi3;


        Eigen::Matrix3d M1;
        Eigen::Matrix3d M2;
        Eigen::Matrix3d M3;

        Eigen::Vector3d h1;
        Eigen::Vector3d h2;
        Eigen::Vector3d h3;


        Eigen::Matrix4d A1 = Eigen::Matrix4d::Zero();
        Eigen::Matrix4d A2 = Eigen::Matrix4d::Zero();
        Eigen::Matrix4d A3 = Eigen::Matrix4d::Zero();

        Eigen::Matrix4d alpha1;
        Eigen::Matrix4d alpha2;
        Eigen::Matrix4d alpha3;

        Eigen::Matrix4d C1;
        Eigen::Matrix4d C2;

        Eigen::Matrix4d Omega;

        Eigen::Vector4d X(0, 0, 0, 1);


        while(t_state.KeepRunning()){
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


                M1 = ::LieAlgebra::skew( Xi1.block<3,1>(0, 0) );
                M2 = ::LieAlgebra::skew( Xi2.block<3,1>(0, 0) );
                M3 = ::LieAlgebra::skew( Xi3.block<3,1>(0, 0) );




                dot_Xi1 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1) * dot_q;
                dot_Xi2 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2) * dot_q;
                dot_Xi3 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3) * dot_q;

                h1 = dot_Xi1.block<3, 1>(0, 0);
                h2 = dot_Xi2.block<3, 1>(0, 0);
                h3 = dot_Xi3.block<3, 1>(0, 0);




                A1.block<3, 4>(0, 0) << M1, h1;
                A2.block<3, 4>(0, 0) << M2, h2;
                A3.block<3, 4>(0, 0) << M3, h3;



                alpha1 = h*A2;
                alpha2 = (A3 - A1)*sqrt(15.0)*h/3.0;
                alpha3 = (A3 - 2*A2 + A1)*10.0*h/3.0;


                C1 = alpha1*alpha2 - alpha2*alpha1;
                C2 = -( alpha1*(2*alpha3 + C1) - (2*alpha3 + C1)*alpha1 )/60.0;

                Omega = alpha1
                        + alpha3/12.0
                        + ( (-20*alpha1 - alpha3 + C1)*(alpha2 + C2)
                            - (alpha2 + C2)*(-20*alpha1 - alpha3 + C1) )/240.0;


                X = Omega.exp() * X;

            }
        }

    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);




    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks("Omega Magnus");



    return 0;
}
