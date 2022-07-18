
#include "benchmark/benchmark.h"

#include "cosserat_rod.hpp"






int main(int argc, char *argv[])
{


//    benchmark::RegisterBenchmark("Forward Kinematics 5pt", [](benchmark::State &t_state){
//        CROSP::CosseratRod rod( 5 );


//        while(t_state.KeepRunning())
//            rod.forwardKinematics();
//    })->Repetitions(10);

//    benchmark::RegisterBenchmark("Forward Kinematics 10pt", [](benchmark::State &t_state){
//        CROSP::CosseratRod rod( 10 );


//        while(t_state.KeepRunning())
//            rod.forwardKinematics();
//    })->Repetitions(10);


//    benchmark::RegisterBenchmark("Forward Kinematics 15pt", [](benchmark::State &t_state){
//        CROSP::CosseratRod rod( 15 );


//        while(t_state.KeepRunning())
//            rod.forwardKinematics();
//    })->Repetitions(10);


    CROSP::CosseratRod rod;
    benchmark::RegisterBenchmark("Forward Kinematics 16pt", [&](benchmark::State &t_state){

        while(t_state.KeepRunning())
            rod.forwardKinematics();
    })->Repetitions(10);


//    benchmark::RegisterBenchmark("Forward Kinematics 20pt", [](benchmark::State &t_state){
//        CROSP::CosseratRod rod( 20 );


//        while(t_state.KeepRunning())
//            rod.forwardKinematics();
//    })->Repetitions(10);

//    benchmark::RegisterBenchmark("Forward Kinematics 25pt", [](benchmark::State &t_state){
//        CROSP::CosseratRod rod( 25 );


//        while(t_state.KeepRunning())
//            rod.forwardKinematics();
//    })->Repetitions(10);

//    benchmark::RegisterBenchmark("Forward Kinematics 30pt", [](benchmark::State &t_state){
//        CROSP::CosseratRod rod( 30 );


//        while(t_state.KeepRunning())
//            rod.forwardKinematics();
//    })->Repetitions(10);

//    benchmark::RegisterBenchmark("Forward Kinematics 35pt", [](benchmark::State &t_state){
//        CROSP::CosseratRod rod( 35 );


//        while(t_state.KeepRunning())
//            rod.forwardKinematics();
//    })->Repetitions(10);


    benchmark::Initialize(&argc, argv);

    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
