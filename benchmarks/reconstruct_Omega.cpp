#include "chebyshev_reconstructor.hpp"

#include <CROSP/CROSP/cosserat_rod.hpp>


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




std::vector<double> getIntermediatePoints()
{
    const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(Nc);

    double begin, end;
    double h, step;

    double x1, x2, x3;

    std::vector<double> intermediate_points;
    const auto number_of_points = 101;
    const auto dX = 1.0/number_of_points;

    for(unsigned int i=0; i<=number_of_points; i++)
        intermediate_points.push_back( i*dX );


    return intermediate_points;


    for(unsigned int i=0; i<Nc-1; i++){


        begin = Chebyshev_points[ (Nc-1) - i];
        end = Chebyshev_points[ (Nc-1) - (i + 1)];

        h = end - begin;
        step = sqrt(15.0)/10.0;

        x1 = begin + (0.5 - step)*h;
        x2 = begin + 0.5*h;
        x3 = begin + (0.5 + step)*h;

        intermediate_points.push_back( x1 );
        intermediate_points.push_back( x2 );
        intermediate_points.push_back( x3 );


    }

    return intermediate_points;
}



void printVector(const std::vector<double> &t_v,
                 const std::string &t_name)
{
    std::cout << t_name << " : ";
    for(const auto &e : t_v)
        std::cout << e << ", ";


    std::cout << "\n\n";
}


int main(int argc, char *argv[])
{

    const ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, ne);

    ::CROSP::CosseratRod rod(polynomial_representation, Nc);

    const auto intermediate_points = getIntermediatePoints();

    ChebyshevReconstructor reconstructor(Nc-2, intermediate_points);



    const Eigen::VectorXd q = Eigen::VectorXd::Random(ne);
    const Eigen::VectorXd dot_q = Eigen::VectorXd::Random(ne);



    //printVector(::Chebyshev::ComputeChebyshevPoints(Nc), "Chebyshev_points");
    printVector(intermediate_points, "intermediate_points");




    rod.updateParameterisation(q, dot_q, 0*dot_q);

    rod.m_idm_integrators->m_angular_velocity->solveSystem();
    const auto Omega_stack = rod.m_idm_integrators->m_angular_velocity->getStackAsMatrix();


    const auto Omega_at_intemediary_points = reconstructor.reconstructVector(Omega_stack.block<3, Nc-2>(0, 1));


    std::cout << "Omega_stack :\n" << Omega_stack.transpose() << "\n\n";

    std::cout << "Omega_at_intemediary_points :\n" << Omega_at_intemediary_points << "\n\n";





    return 0;
}
