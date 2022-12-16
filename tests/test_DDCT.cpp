#include "chebyshev_reconstructor.hpp"

#include <vector>
#include <algorithm>

#include <iostream>

#include <boost/math/special_functions/chebyshev.hpp>

const unsigned int Nc = 8;

void printVector(const std::vector<double> &t_v,
                 const std::string &t_name)
{
    std::cout << t_name << " : ";
    for(const auto &e : t_v)
        std::cout << e << ", ";


    std::cout << "\n\n";
}



double Tn_xm(const unsigned int n,
             const unsigned int m)
{
    return cos(M_PI*n*m/(Nc-1));
}


double am(const unsigned int t_m,
          const std::vector<double> u_stack,
          const std::vector<double> Chebyshev_points)
{
    const auto pm = [](const unsigned int m){
        if(m == 0 || m == Nc-1)
            return 1.0;

        return 2.0;
    };

    double range = 0.5 * (pow(-1, t_m) * u_stack[0] + u_stack[Nc-1]);

    double sum = 0;
    for(unsigned int n=1; n<Nc; n++){
        sum += u_stack[n] * Tn_xm(t_m, n);
    }

    double a_m = pm(t_m)/(Nc-1)*( range + sum );

    return a_m;

}


int main(int argc, char *argv[])
{


    std::vector<double> Chebyshev_points(Nc);
    for(unsigned int j=0; j<Nc; j++){
        auto xj = cos( M_PI * j / (Nc-1.0) );
        std::cout << "xj : " << xj << "\n";

        Chebyshev_points[j] = 0.5*(1 + xj);
    }


    std::vector<double> f_constant(Nc);
    std::vector<double> f_x(Nc);

    for(unsigned int j=0; j<Nc; j++){
        f_constant[j] = 3;
        f_x[j] = Chebyshev_points[(Nc-1)-j];
    }




    const unsigned int number_of_points = 21;
    std::vector<double> interpolation_points(number_of_points);



    const double step = 1.0/(number_of_points-1);
    for(unsigned int j=0; j<number_of_points; j++)
        interpolation_points[j] = j*step;





    printVector(Chebyshev_points, "Chebyshev_points");

    printVector(interpolation_points, "interpolation_points");


    std::vector<double> u_stack(number_of_points);

    for(unsigned int n=0; n<number_of_points; n++){
        u_stack[n] = 0;
        for(unsigned int m=0; m<Nc; m++){
            u_stack[n] += am(m, f_constant, Chebyshev_points) * Tn_xm(m, n);
        }
    }

    printVector(u_stack, "u_stack");




//    ChebyshevReconstructor reconstructor(Nc, interpolation_points, 1);

//    std::cout << "DDCT : \n" << (Nc/2.0)*reconstructor.m_DDCT << "\n\n";


//    auto res_constant = reconstructor.reconstructVector( f_constant );
//    auto res_x = reconstructor.reconstructVector( f_x );


//    std::cout << "res_constant : \n" << res_constant.transpose() << "\n\n";

//    std::cout << "res_x : \n" << res_x.transpose() << "\n\n";



    return 0;
}
