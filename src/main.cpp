#include "CROSP/CROSP/cosserat_rod.hpp"


int main(int argc, char *argv[])
{

    ::CROSP::polynomial_representation::PolynomialRepresentation poly(5);

    ::CROSP::rod_properties::MaterialProperties material(200, 100, 10);

    ::CROSP::CosseratRod rod(poly, material);

    ::CROSP::CosseratRod rod2(15);

    ::CROSP::CosseratRod rod3;
    return 0;
}





