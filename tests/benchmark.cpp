
#include "benchmark/benchmark.h"

#include "CROSP/cosserat_rod.hpp"



int main(int argc, char *argv[])
{


    benchmark::RegisterBenchmark("Forward Kinematics 10", [](benchmark::State &t_state){

        CROSP::CosseratRod rod(10);

        while(t_state.KeepRunning())
            rod.forwardKinematics();
    })->Repetitions(10);


    benchmark::RegisterBenchmark("Forward Kinematics 15", [](benchmark::State &t_state){

        CROSP::CosseratRod rod(15);

        while(t_state.KeepRunning())
            rod.forwardKinematics();
    })->Repetitions(10);

    benchmark::RegisterBenchmark("Forward Kinematics 20", [](benchmark::State &t_state){

        CROSP::CosseratRod rod(20);

        while(t_state.KeepRunning())
            rod.forwardKinematics();
    })->Repetitions(10);

    benchmark::RegisterBenchmark("Forward Kinematics 25", [](benchmark::State &t_state){

        CROSP::CosseratRod rod(25);

        while(t_state.KeepRunning())
            rod.forwardKinematics();
    })->Repetitions(10);

    benchmark::RegisterBenchmark("Forward Kinematics 30", [](benchmark::State &t_state){

        CROSP::CosseratRod rod(30);

        while(t_state.KeepRunning())
            rod.forwardKinematics();
    })->Repetitions(10);

    benchmark::RegisterBenchmark("Forward Kinematics 35", [](benchmark::State &t_state){

        CROSP::CosseratRod rod(35);

        while(t_state.KeepRunning())
            rod.forwardKinematics();
    })->Repetitions(10);


    benchmark::Initialize(&argc, argv);

    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
