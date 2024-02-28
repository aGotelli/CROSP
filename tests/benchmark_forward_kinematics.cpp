
#include "benchmark/benchmark.h"

#include "CROSP/CROSP/cosserat_rod.hpp"



int main(int argc, char *argv[])
{



    std::vector<unsigned int> Chebyshev_points_staks {
        10,
        15,
        20,
        25,
        30,
        35,
        40
    };


    for(const auto number_of_Chebyshev_points : Chebyshev_points_staks) {
        benchmark::RegisterBenchmark("Forward Kinematics with Nc =", [](benchmark::State &t_state){



            auto strain_param =
                    std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(t_state.range(0));

            CROSP::CosseratRod rod(strain_param);


            while(t_state.KeepRunning())
                rod.forwardKinematics();
        })->Repetitions(10)->Arg(number_of_Chebyshev_points);
    }



    benchmark::Initialize(&argc, argv);

    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
