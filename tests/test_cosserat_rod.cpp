#include "CROSP/CROSP/cosserat_rod.hpp"

#include "utilities/Eigen/eigen_io.hpp"





int main()
{
    constexpr unsigned int number_of_Chebyshev_points = 17;
    CROSP::CosseratRod rod;


    writeToFile("Kee_cpp", rod.m_rod_properties->m_Kee, "/home/andrea/Desktop/PhD/PhD_development/strain_approach/MATLAB/Dyn_Essai_release_Beam_Andrea/data_from_cpp");

    return 0;
}
