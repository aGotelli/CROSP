#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "CROSP/static_initialiser/static_initialiser.hpp"

#include "utilities/Eigen/eigen_io.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/range/algorithm.hpp>

#include <array>



constexpr std::array<bool, 6> admitted_deformations = {
    false,
    true,
    false,

    false,
    false,
    false
};

constexpr unsigned int ne = 4;





std::vector<Eigen::Matrix3d> solutionSpectral(::CROSP::CosseratRod t_rod,
                                 const Eigen::VectorXd &t_q)
{
    const double number_of_Chebyshev_points = t_rod.m_strain_parameterisation->m_number_of_Chebyshev_points;

    t_rod.updateParameterisation(t_q, 0*t_q, 0*t_q);

    t_rod.m_idm_integrators->m_quaternion->solveSystem();


    std::vector<Eigen::Matrix3d> R_stack;

    for(unsigned int i=0; i<number_of_Chebyshev_points; i++){

        const auto Q = t_rod.m_idm_integrators->m_quaternion->getStateAtPoint(i);

        R_stack.push_back( Eigen::Quaterniond(Q[0],
                           Q[1],
                           Q[2],
                           Q[3]).toRotationMatrix() );
    }


    return R_stack;

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


std::vector<Eigen::Matrix3d> solutionODE(const ::CROSP::CosseratRod &rod,
                                         const Eigen::VectorXd &q,
                                         std::vector<double> t_Chebyshev_points)
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
    const double dX = 0.001;

    std::vector<Eigen::Matrix3d> R_stack;

    const auto points_reversed = boost::range::reverse(t_Chebyshev_points);

    boost::numeric::odeint::integrate_times(stepper(),
                                            [&](const state_type& t_y, state_type& t_dyds, const double t_X){

        Xi = rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(t_X) * q;

        k = Xi.block<3,1>(0,0);

        t_dyds = 0.5*getA(k)*t_y;

    }, y, points_reversed.begin(), points_reversed.end(), dX, [&](const state_type& t_y, const double){



        R_stack.push_back( Eigen::Quaterniond(t_y[0],
                           t_y[1],
                           t_y[2],
                           t_y[3]).toRotationMatrix() );
    });



    return R_stack;

}



std::vector<Eigen::Matrix3d> solutionMagnusOptimised(const ::CROSP::CosseratRod &rod,
                                                     const Eigen::VectorXd &q)
{

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


//        std::cout << "at i = " << i << "\n";
//        std::cout << "[x1, x2, x3] : [" << x1 << ", " << x2 << ", " << x3 << "]\n" << std::endl;


        B_Phi_stack[i] = {rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x1),
                          rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x2),
                          rod.m_polynomial_representation.m_B * rod.m_polynomial_representation.getPhi(x3)};

        h_stack[i] = h;

    }



    ::LieAlgebra::Vector6d Xi1;
    ::LieAlgebra::Vector6d Xi2;
    ::LieAlgebra::Vector6d Xi3;




    std::vector<QuadrarureObserver> K_at_quadrature_points(number_of_Chebyshev_points-1);
    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){

        Xi1 = B_Phi_stack[i].v1 * q;
        Xi2 = B_Phi_stack[i].v2 * q;
        Xi3 = B_Phi_stack[i].v3 * q;


        K_at_quadrature_points[i] = { Xi1.block<3,1>(0, 0),
                                      Xi2.block<3,1>(0, 0),
                                      Xi3.block<3,1>(0, 0) };
    }

    Eigen::Matrix3d A1;
    Eigen::Matrix3d A2;
    Eigen::Matrix3d A3;

    Eigen::Matrix3d alpha1;
    Eigen::Matrix3d alpha2;
    Eigen::Matrix3d alpha3;

    Eigen::Matrix3d C1;
    Eigen::Matrix3d C2;

    Eigen::Matrix3d Omega;

    std::vector<Eigen::Matrix3d> exp_Omega(number_of_Chebyshev_points-1);




    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++){


        A1 = ::LieAlgebra::skew( K_at_quadrature_points[i].v1 );
        A2 = ::LieAlgebra::skew( K_at_quadrature_points[i].v2 );
        A3 = ::LieAlgebra::skew( K_at_quadrature_points[i].v3 );


        alpha1 = h_stack[i]*A2;
        alpha2 = (A3 - A1)*sqrt(15.0)*h_stack[i]/3.0;
        alpha3 = (A3 - 2*A2 + A1)*10.0*h_stack[i]/3.0;


        C1 = alpha1*alpha2 - alpha2*alpha1;
        C2 = -( alpha1*(2*alpha3 + C1) - (2*alpha3 + C1)*alpha1 )/60.0;

        Omega = alpha1
                + alpha3/12.0
                + ( (-20*alpha1 - alpha3 + C1)*(alpha2 + C2)
                    - (alpha2 + C2)*(-20*alpha1 - alpha3 + C1) )/240.0;

        exp_Omega[i] = ::LieAlgebra::expRodigues(Omega);




    }


    std::vector<Eigen::Matrix3d> R(number_of_Chebyshev_points);
    R[0] = Eigen::Matrix3d::Identity();

    for(unsigned int i=0; i<number_of_Chebyshev_points-1; i++)
        R[i+1] = R[i] * exp_Omega[i];


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




    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);



    std::vector<unsigned int> Nc_stack {
        5,
        11,
        15,
//        21,
//        25,
//        31
    };


    Eigen::VectorXd q = Eigen::VectorXd::Zero(ne);
    q.setRandom();




    Eigen::MatrixXd errors_along_rod_OSNI(Nc_stack.size(), Nc_stack[Nc_stack.size()-1]);
    Eigen::MatrixXd errors_along_rod_Magnus(Nc_stack.size(), Nc_stack[Nc_stack.size()-1]);

    for(unsigned int j=0;const auto Nc : Nc_stack){
        ::CROSP::CosseratRod rod(polynomial_representation, Nc);

        const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(Nc);

        const auto R_stack_ODE              = solutionODE(rod, q, Chebyshev_points);
        const auto R_stack_OSNI             = solutionSpectral(rod, q);
        const auto R_stack_Magnus_optimised = solutionMagnusOptimised(rod, q);
        //const auto R_Cayley           = solutionCayley(rod, q);



        for(unsigned int i=0; i<Nc; i++){
//            std::cout << "\n\n" << "iter : " << i << "\n";
//            std::cout << "R ODE : \n" << R_stack_ODE[i] << "\n";
//            std::cout << "R OSNI : \n" << R_stack_OSNI[(Nc - 1) - i] << "\n";
//            std::cout << std::endl;
            auto error = ::LieAlgebra::logSO3(R_stack_ODE[i], R_stack_OSNI[(Nc - 1) - i]);
            errors_along_rod_OSNI(j, i) = error[1];


            error = ::LieAlgebra::logSO3(R_stack_ODE[i], R_stack_Magnus_optimised[i]);
            errors_along_rod_Magnus(j, i) = error[1];
        }

        j++;

    }

//    std::cout << "errors_along_rod_OSNI : \n" << errors_along_rod_OSNI << "\n\n";
//    std::cout << "errors_along_rod_Magnus : \n" << errors_along_rod_Magnus << "\n\n";

    writeToFile("errors_along_rod_OSNI", errors_along_rod_OSNI);
    writeToFile("errors_along_rod_Magnus", errors_along_rod_Magnus);




}





int main(int argc, char *argv[])
{

    compareSolutions();



    return 0;
}
