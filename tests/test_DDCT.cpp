#include "chebyshev_reconstructor.hpp"

#include <vector>
#include <algorithm>

#include <iostream>

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
    const unsigned int Nc = 8;

    std::vector<double> Chebyshev_points(Nc);
    for(unsigned int j=0; j<Nc; j++){
        auto xj = cos(M_PI*(2*j + 1)/(2*Nc));

        Chebyshev_points[j] = 0.5*(1 + xj);
    }


    Eigen::MatrixXd f_constant(1, Nc);
    Eigen::MatrixXd f_x(1, Nc);

    for(unsigned int j=0; j<Nc; j++){
        f_constant(0, j) = 3;
        f_x(0, j) = Chebyshev_points[(Nc-1)-j];
    }


    std::vector<double> interpolation_points;

//    for(unsigned int j=0; j<Nc-1; j++){

//        auto begin = Chebyshev_points[ (Nc-1) - j];
//        auto end = Chebyshev_points[ (Nc-1) - (j + 1)];

//        auto h = end - begin;
//        auto step = sqrt(15.0)/10.0;

//        auto x1 = begin + (0.5 - step)*h;
//        auto x2 = begin + 0.5*h;
//        auto x3 = begin + (0.5 + step)*h;

//        interpolation_points.push_back( x1 );
//        interpolation_points.push_back( x2 );
//        interpolation_points.push_back( x3 );
//    }
    const unsigned int number_of_points = 21;
    const double step = 1.0/(number_of_points-1);
    for(unsigned int j=0; j<number_of_points; j++){


        interpolation_points.push_back( j*step );

    }




    printVector(Chebyshev_points, "Chebyshev_points");

    printVector(interpolation_points, "interpolation_points");


    ChebyshevReconstructor reconstructor(Nc, interpolation_points, 1);

    std::cout << "DDCT : \n" << (Nc/2.0)*reconstructor.m_DDCT << "\n\n";


    auto res_constant = reconstructor.reconstructVector( f_constant );
    auto res_x = reconstructor.reconstructVector( f_x );


    std::cout << "res_constant : \n" << res_constant.transpose() << "\n\n";

    std::cout << "res_x : \n" << res_x.transpose() << "\n\n";



    return 0;
}
