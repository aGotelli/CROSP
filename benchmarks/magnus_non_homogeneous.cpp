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

constexpr unsigned int Nc = 31;



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
        std::cout << "i = " << i <<"\n";

        begin = Chebyshev_points[ (number_of_Chebyshev_points-1) - i];
        end = Chebyshev_points[ (number_of_Chebyshev_points-1) - (i + 1)];


//        std::cout << "begin : " << begin << "\n";
//        std::cout << "end : " << end << "\n";

        h = end - begin;
        step = sqrt(15.0)/10.0;

        x1 = begin + (0.5 - step)*h;
        x2 = begin + 0.5*h;
        x3 = begin + (0.5 + step)*h;

        const Eigen::MatrixXd BPhi1 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1);
        const Eigen::MatrixXd BPhi2 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2);
        const Eigen::MatrixXd BPhi3 = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3);

//        std::cout << "[x1, x2, x3] : [" << x1 << ", " << x2 << ", " << x3 << "]\n" << std::endl;

//        std::cout << "BPhi(x1) : \n" << BPhi1 << "\n\n";
//        std::cout << "BPhi(x2) : \n" << BPhi2 << "\n\n";
//        std::cout << "BPhi(x3) : \n" << BPhi3 << "\n\n";
        std::cout << "\n\n\n" << std::endl;




        Xi1 = BPhi1 * t_q;
        Xi2 = BPhi2 * t_q;
        Xi3 = BPhi3 * t_q;


        M1 = ::LieAlgebra::skew( Xi1.block<3,1>(0, 0) );
        M2 = ::LieAlgebra::skew( Xi2.block<3,1>(0, 0) );
        M3 = ::LieAlgebra::skew( Xi3.block<3,1>(0, 0) );


//        std::cout << "q : " << t_q.transpose() << "\n" << std::endl;

//        std::cout << "K1 : \n" << Xi1.block<3,1>(0, 0) << "\n";

//        std::cout << "M1 : \n" << M1 << "\n";
//        std::cout << "M2 : \n" << M2 << "\n";
//        std::cout << "M3 : \n" << M3 << "\n";
//        std::cout << "\n\n" << std::endl;




        dot_Xi1 = BPhi1 * t_dot_q;
        dot_Xi2 = BPhi2 * t_dot_q;
        dot_Xi3 = BPhi3 * t_dot_q;

        h1 = dot_Xi1.block<3, 1>(0, 0);
        h2 = dot_Xi2.block<3, 1>(0, 0);
        h3 = dot_Xi3.block<3, 1>(0, 0);

        std::cout << "h1 : \n" << h1 << "\n";
        std::cout << "h2 : \n" << h2 << "\n";
        std::cout << "h3 : \n" << h3 << "\n";
        std::cout << "\n\n" << std::endl;




        A1.block<3, 4>(0, 0) << M1, h1;
        A2.block<3, 4>(0, 0) << M2, h2;
        A3.block<3, 4>(0, 0) << M3, h3;



        alpha1 = h*A2;
        alpha2 = (A3 - A1)*sqrt(15.0)*h/3.0;
        alpha3 = (A3 - 2*A2 + A1)*10.0*h/3.0;

//        std::cout << "alpha1 : \n" << alpha1 << "\n";
//        std::cout << "alpha2 : \n" << alpha2 << "\n";
//        std::cout << "alpha3 : \n" << alpha3 << "\n";
//        std::cout << "\n\n" << std::endl;


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




Eigen::MatrixXd solutionMagnusOptimised(const ::CROSP::CosseratRod &rod,
                                        const Eigen::VectorXd &q,
                                        const Eigen::VectorXd &dot_q)
{

    std::cout << "\n\n\n\n\n\n" << "Magnus optimised" << "\n\n\n\n\n\n";
    const double number_of_Chebyshev_points = rod.m_strain_parameterisation->m_number_of_Chebyshev_points;
    const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(number_of_Chebyshev_points);


    double x1, x2, x3;

    struct QuadrarureObserver{
        Eigen::MatrixXd v1;
        Eigen::MatrixXd v2;
        Eigen::MatrixXd v3;
    };

    std::vector<double> h_stack(number_of_Chebyshev_points-1);
    std::vector<QuadrarureObserver> B_Phi_stack(number_of_Chebyshev_points-1);
    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){
//        std::cout << "at i = " << i << "\n";


        const auto begin = Chebyshev_points[ (number_of_Chebyshev_points-1) - i];
        const auto end = Chebyshev_points[ (number_of_Chebyshev_points-1) - (i + 1)];


//        std::cout << "begin : " << begin << "\n";
//        std::cout << "end : " << end << "\n";

        const auto h = end - begin;
        const auto step = sqrt(15.0)/10.0;

        x1 = begin + (0.5 - step)*h;
        x2 = begin + 0.5*h;
        x3 = begin + (0.5 + step)*h;

//        std::cout << "[x1, x2, x3] : [" << x1 << ", " << x2 << ", " << x3 << "]\n" << std::endl;




        B_Phi_stack[i] = {rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1),
                          rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2),
                          rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3)};

        h_stack[i] = h;

    }



    ::LieAlgebra::Vector6d Xi1;
    ::LieAlgebra::Vector6d Xi2;
    ::LieAlgebra::Vector6d Xi3;


    ::LieAlgebra::Vector6d dot_Xi1;
    ::LieAlgebra::Vector6d dot_Xi2;
    ::LieAlgebra::Vector6d dot_Xi3;




    std::vector<QuadrarureObserver> K_at_quadrature_points(number_of_Chebyshev_points-1);
    std::vector<QuadrarureObserver> dot_K_at_quadrature_points(number_of_Chebyshev_points-1);
    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){
//        std::cout << "at i = " << i << "\n";

//        std::cout << "BPhi(x1) : \n" << B_Phi_stack[i].v1 << "\n\n";
//        std::cout << "BPhi(x2) : \n" << B_Phi_stack[i].v2 << "\n\n";
//        std::cout << "BPhi(x3) : \n" << B_Phi_stack[i].v3 << "\n\n";
//        std::cout << "\n\n\n" << std::endl;

//        std::cout << "q : " << q.transpose() << "\n" << std::endl;

        Xi1 = B_Phi_stack[i].v1 * q;
        Xi2 = B_Phi_stack[i].v2 * q;
        Xi3 = B_Phi_stack[i].v3 * q;




        dot_Xi1 = B_Phi_stack[i].v1 * dot_q;
        dot_Xi2 = B_Phi_stack[i].v2 * dot_q;
        dot_Xi3 = B_Phi_stack[i].v3 * dot_q;


//        std::cout << "xi 1 : \n" << Xi1 << "\n" << std::endl;

        K_at_quadrature_points[i] = { Xi1.block<3,1>(0, 0),
                                      Xi2.block<3,1>(0, 0),
                                      Xi3.block<3,1>(0, 0) };

        dot_K_at_quadrature_points[i] = { dot_Xi1.block<3,1>(0, 0),
                                          dot_Xi2.block<3,1>(0, 0),
                                          dot_Xi3.block<3,1>(0, 0) };
    }

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

    std::vector<Eigen::Matrix4d> Omega(number_of_Chebyshev_points-1);




    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){

        std::cout << "i = " << i <<"\n";


        M1 = ::LieAlgebra::skew( K_at_quadrature_points[i].v1 );
        M2 = ::LieAlgebra::skew( K_at_quadrature_points[i].v2 );
        M3 = ::LieAlgebra::skew( K_at_quadrature_points[i].v3 );

//        std::cout << "K1 : \n" << K_at_quadrature_points[i].v1 << "\n";

//        std::cout << "M1 : \n" << M1 << "\n";
//        std::cout << "M2 : \n" << M2 << "\n";
//        std::cout << "M3 : \n" << M3 << "\n";
//        std::cout << "\n\n" << std::endl;


        h1 = dot_K_at_quadrature_points[i].v1;
        h2 = dot_K_at_quadrature_points[i].v2;
        h3 = dot_K_at_quadrature_points[i].v3;

        std::cout << "h1 : \n" << h1 << "\n";
        std::cout << "h2 : \n" << h2 << "\n";
        std::cout << "h3 : \n" << h3 << "\n";
        std::cout << "\n\n" << std::endl;


        A1.block<3, 4>(0, 0) << M1, h1;
        A2.block<3, 4>(0, 0) << M2, h2;
        A3.block<3, 4>(0, 0) << M3, h3;


        alpha1 = h_stack[i]*A2;
        alpha2 = (A3 - A1)*sqrt(15.0)*h_stack[i]/3.0;
        alpha3 = (A3 - 2*A2 + A1)*10.0*h_stack[i]/3.0;

//        std::cout << "alpha1 : \n" << alpha1 << "\n";
//        std::cout << "alpha2 : \n" << alpha2 << "\n";
//        std::cout << "alpha3 : \n" << alpha3 << "\n";
//        std::cout << "\n\n" << std::endl;


        C1 = alpha1*alpha2 - alpha2*alpha1;
        C2 = -( alpha1*(2*alpha3 + C1) - (2*alpha3 + C1)*alpha1 )/60.0;

        Omega[i] = alpha1
                    + alpha3/12.0
                    + ( (-20*alpha1 - alpha3 + C1)*(alpha2 + C2)
                        - (alpha2 + C2)*(-20*alpha1 - alpha3 + C1) )/240.0;

    }


    std::vector<Eigen::Vector4d> X_stack(number_of_Chebyshev_points);
    X_stack[0] << 0, 0, 0, 1;

    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++)
        X_stack[i+1] = Omega[i].exp() * X_stack[i];


    return X_stack[number_of_Chebyshev_points-1].block<3, 1>(0, 0);

}






void compareSolutions()
{


    const ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation, Nc);




    const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);
    const Eigen::VectorXd dot_q = Eigen::VectorXd::Random(ne);


    const auto Omega_X1_OSNI   = solutionSpectral(rod, q, dot_q);
    const auto Omega_X1_Magnus = solutionMagnus(rod, q, dot_q);
    const auto Omega_X1_Magnus_optimised = solutionMagnusOptimised(rod, q, dot_q);

    std::cout << "With \n" <<
                 "      q = [" << q.transpose() << "]\n" <<
                 "  dot_q = [" << dot_q.transpose() << "]\n\n";


    std::cout << "Solutions :\n" <<
                 "  OSNI\n" <<
                 "      Omega : \n" << Omega_X1_OSNI << "\n\n" <<
                 "  Magnus\n" <<
                 "      Omega : \n" << Omega_X1_Magnus << "\n\n" <<
                 "  Magnus optimised\n" <<
                 "      Omega : \n" << Omega_X1_Magnus_optimised << "\n\n";

    std::cout << "\n\n\n";

    const auto error_OSNI_Magnus = Omega_X1_Magnus - Omega_X1_OSNI;
    const auto error_OSNI_Magnus_optimised = Omega_X1_Magnus_optimised - Omega_X1_OSNI;

    std::cout << "Errors :\n" <<
                 "  OSNI-Magnus\n" <<
                 error_OSNI_Magnus << "\n\n" <<
                 "  OSNI-Magnus optimised\n" <<
                 error_OSNI_Magnus_optimised << "\n\n";


}





int main(int argc, char *argv[])
{


    srand(0);

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



/*
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
*/




    ::benchmark::RegisterBenchmark("Omega Magnus optimised", [&](::benchmark::State &t_state) {

        const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);
        const Eigen::VectorXd dot_q = Eigen::VectorXd::Random(ne);

        ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

        ::CROSP::CosseratRod rod(polynomial_representation, Nc);

        rod.updateParameterisation(q, dot_q, 0*dot_q);


        const double number_of_Chebyshev_points = rod.m_strain_parameterisation->m_number_of_Chebyshev_points;
        const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(number_of_Chebyshev_points);


        double x1, x2, x3;

        struct QuadrarureObserver{
            Eigen::MatrixXd v1;
            Eigen::MatrixXd v2;
            Eigen::MatrixXd v3;
        };

        std::vector<double> h_stack(number_of_Chebyshev_points-1);
        std::vector<QuadrarureObserver> B_Phi_stack(number_of_Chebyshev_points-1);
        for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){


            const auto begin = Chebyshev_points[ (number_of_Chebyshev_points-1) - i];
            const auto end = Chebyshev_points[ (number_of_Chebyshev_points-1) - (i + 1)];


            const auto h = end - begin;
            const auto step = sqrt(15.0)/10.0;

            x1 = begin + (0.5 - step)*h;
            x2 = begin + 0.5*h;
            x3 = begin + (0.5 + step)*h;


            B_Phi_stack[i] = {rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1),
                              rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2),
                              rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3)};

            h_stack[i] = h;

        }



        ::LieAlgebra::Vector6d Xi1;
        ::LieAlgebra::Vector6d Xi2;
        ::LieAlgebra::Vector6d Xi3;


        ::LieAlgebra::Vector6d dot_Xi1;
        ::LieAlgebra::Vector6d dot_Xi2;
        ::LieAlgebra::Vector6d dot_Xi3;




        std::vector<QuadrarureObserver> K_at_quadrature_points(number_of_Chebyshev_points-1);
        std::vector<QuadrarureObserver> dot_K_at_quadrature_points(number_of_Chebyshev_points-1);
        for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){

            Xi1 = B_Phi_stack[i].v1 * q;
            Xi2 = B_Phi_stack[i].v2 * q;
            Xi3 = B_Phi_stack[i].v3 * q;




            dot_Xi1 = B_Phi_stack[i].v1 * dot_q;
            dot_Xi2 = B_Phi_stack[i].v2 * dot_q;
            dot_Xi3 = B_Phi_stack[i].v3 * dot_q;


            K_at_quadrature_points[i] = { Xi1.block<3,1>(0, 0),
                                          Xi2.block<3,1>(0, 0),
                                          Xi3.block<3,1>(0, 0) };

            dot_K_at_quadrature_points[i] = { dot_Xi1.block<3,1>(0, 0),
                                              dot_Xi2.block<3,1>(0, 0),
                                              dot_Xi3.block<3,1>(0, 0) };
        }

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

        std::vector<Eigen::Matrix4d> Omega(number_of_Chebyshev_points-1);


        std::vector<Eigen::Vector4d> X_stack(number_of_Chebyshev_points);
        X_stack[0] << 0, 0, 0, 1;


        while(t_state.KeepRunning()){
            for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){


                M1 = ::LieAlgebra::skew( K_at_quadrature_points[i].v1 );
                M2 = ::LieAlgebra::skew( K_at_quadrature_points[i].v2 );
                M3 = ::LieAlgebra::skew( K_at_quadrature_points[i].v3 );


                h1 = dot_K_at_quadrature_points[i].v1;
                h2 = dot_K_at_quadrature_points[i].v2;
                h3 = dot_K_at_quadrature_points[i].v3;


                A1.block<3, 4>(0, 0) << M1, h1;
                A2.block<3, 4>(0, 0) << M2, h2;
                A3.block<3, 4>(0, 0) << M3, h3;


                alpha1 = h_stack[i]*A2;
                alpha2 = (A3 - A1)*sqrt(15.0)*h_stack[i]/3.0;
                alpha3 = (A3 - 2*A2 + A1)*10.0*h_stack[i]/3.0;


                C1 = alpha1*alpha2 - alpha2*alpha1;
                C2 = -( alpha1*(2*alpha3 + C1) - (2*alpha3 + C1)*alpha1 )/60.0;

                Omega[i] = alpha1
                            + alpha3/12.0
                            + ( (-20*alpha1 - alpha3 + C1)*(alpha2 + C2)
                                - (alpha2 + C2)*(-20*alpha1 - alpha3 + C1) )/240.0;

            }



            for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++)
                X_stack[i+1] = Omega[i].exp() * X_stack[i];
        }

    })->Repetitions(repetitions)->Unit(::benchmark::kMicrosecond);




    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
