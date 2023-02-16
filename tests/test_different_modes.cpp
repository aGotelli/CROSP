#include "CROSP/CROSP/cosserat_rod.hpp"










int main(int argc, char *argv[])
{

    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation;

    std::cout << "Case with standand default constructor : \n"
                 "  Phi(X=0) : \n" <<
                 polynomial_representation.getPhi(0) << "\n\n";



    std::array<bool, 6> admitted_deformations{
        true,
        true,
        true,
        false,
        false,
        false
    };

    unsigned int ne = 3;


    ::CROSP::polynomial_representation::PolynomialRepresentation poly_with_admitted_def(admitted_deformations);

    std::cout << "Case admitted deformations : \n"
                 "  Phi(X=0) : \n" <<
                 poly_with_admitted_def.getPhi(0) << "\n\n";



    admitted_deformations = {
            true,
            true,
            true,
            true,
            false,
            false
        };

    std::vector<unsigned int> number_of_modes_stack {
        2,
        5,
        5,
        1
    };


    ::CROSP::polynomial_representation::PolynomialRepresentation poly_with_admitted_def_and_modes(admitted_deformations, number_of_modes_stack);

    std::cout << "Case admitted deformations : \n"
                 "  Phi(X=0) : \n" <<
                 poly_with_admitted_def_and_modes.getPhi(0) << "\n\n";







    //  Now use it in rod properties
    ::CROSP::rod_properties::RodProperties rod_properties(poly_with_admitted_def_and_modes);

    return 0;
}
