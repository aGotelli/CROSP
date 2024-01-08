clear all
clc

Theta =[ 
-0.85911;
-2.44664;
 0.51629];

Delta_Theta = [
  0.59688;
 0.823295;
-0.604897];


Delta_T_SO3_cpp = [
 0.455097 -0.250313 0.0529811;
   0.3919   0.63185 -0.333649;
-0.0828588   0.248561   0.439216;
];



Delta_TSO3 = Delta_T_SO3_(Theta,Delta_Theta);

error = Delta_T_SO3_cpp - Delta_TSO3