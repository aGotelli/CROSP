clear all
clc

addpath('./outils')
addpath('./IDM')
addpath('./methode spectral')
Config.option=odeset('RelTol',10^(-8),'AbsTol',10^(-8));

% Variable pour visualisation du transitoire
% 0 pas de visualisation
% 1 visualisation
Config.Trace_transitoire = 0;
    
% Initialisation des variables pour le tracé
r_dot_trace_save = [];
T_trace_save = [];
x_trace_save = [];
y_trace_save = [];
z_trace_save = [];
           
% Initialisation boucle de temps
Config.dt    = 0.005;
Config.t_fin = 10.0;
Config.time  = 0:Config.dt:Config.t_fin;
Config.Nt    = length(Config.time);


% longueur de la poutre
Config.L = 2.41;
% Pas de discrétisation du corps
Config.dX = Config.L/101;
% Construction du corps et de ses paramètres
Config.X  = 0:Config.dX:Config.L;
% Nombre de noeuds discrétisant le corps
Config.Nx = length(Config.X);

Config.V_a = load("../test_different_modes/deformations_stack.csv");
% Définition de la taille de la base pour chaque variable
Const.dim_base_k = load("../test_different_modes/ne_stack.csv");
% Calcul de la taille de q 
Const.dim_base   = Config.V_a'*Const.dim_base_k;


% Valeur de la gravité
Const.Gamma_g = 9.81;


% Parametre constants de la simulation
% Valeurs physique de la poutre 
adapted_parametre_constant_Compare_Dominique_F_ABS;



Const.B = load("../test_different_modes/B.csv");
Const.B_bar = load("../test_different_modes/B_bar.csv");

% Initialisation du temps à 0
time = 0;

% Position et orientation initiales de la tete
r_0 = [0;0;0];
q_0 = [1 0 0 0]';
%q_0 = [0.7071068 0 0.7071068 0]';
eta_0 = [0 0 0 0 0 0]'; % vitesse de la tête
eta_dot_0 = [0 0 0 0 0 0]';

% ------------ Calcul des gain elastique et ammortissement --------------- %

% ------ Initialisation ----- %
Kee_0 = zeros(Const.dim_base*Const.dim_base,1);
Dee_0 = zeros(Const.dim_base*Const.dim_base,1);

% ------- Calcul ------ %
CI = [Kee_0; Dee_0];
[X Y1] = ode45(@(t,y) elastic_gain(t,y,time,Const,Config),Config.X,CI,Config.option); % Solve ODE

% ------- Affectation des résultats -----------%
Kee_vec = Y1(end,1:Const.dim_base*Const.dim_base);
Dee_vec = Y1(end,1+Const.dim_base*Const.dim_base:2*Const.dim_base*Const.dim_base);

Kee = reshape(Kee_vec, [Const.dim_base Const.dim_base]);
Dee = reshape(Dee_vec, [Const.dim_base Const.dim_base]);

Const.Kee = Kee;
% Const.Dee = Dee;
Const.Dee = Const.mu*Kee;

% ------------------------------------------------------- %

% Initialisation de la base q
q0 = zeros(Const.dim_base,1);

% Définiton de la variation de la force en bout 
%C_1 = [0:0.1:0.5];



% Const.q         = Y2_s(end,:)';
% Const.q         =[  -0.2; 0.2; -0.2];
% Const.q_dot     = zeros(Const.dim_base,1);
% Const.q_dot_dot = zeros(Const.dim_base,1);

Const.q = load("../test_different_modes/q.csv");
Const.q_dot = load("../test_different_modes/dot_q.csv");
Const.q_dot_dot = load("../test_different_modes/ddot_q.csv");

display(Const.q)
display(Const.q_dot)
display(Const.q_dot_dot)

Kee_cpp = load("../test_different_modes/Kee.csv");
Dee_cpp = load("../test_different_modes/Dee.csv");

error_Kee = Kee_cpp - Kee;
error_Kee_norm = norm(error_Kee);
display(Kee_cpp)
display(Kee)
display(error_Kee_norm)





error_Dee = Dee_cpp - Dee;
error_Dee_norm = norm(error_Dee);
display(error_Dee_norm)



Const.eta_dot   = zeros(6,1);
Const.eta_dot_0   = zeros(6,1);

Const.r_0 = r_0;
Const.q_0 = q_0; 
Const.eta_0 = eta_0;
Const.F1 = zeros(6,1);
CI = [Const.q;Const.q_dot];



a = 400;
b = 160000;
% a = 0;
% b = 0;

time = 0;
qn_0 = q_0;
rn_0 = r_0;


[QX, rX, eta_X, eta_dot_X, Lambda_X, Qa_X, ...
    Delta_zeta_X, Delta_eta_X, Delta_eta_dot_X, Delta_Lambda_X, Delta_Qa_X, J] = ...
        TIDM_spectral_output_all_data(qn_0,rn_0,q_0,r_0,eta_0,eta_dot_0,a,b,time,Const,Config);


Q_cpp = load("../test_different_modes/Q_stack.csv");
r_cpp = load("../test_different_modes/r_stack.csv");
Omega_cpp = load("../test_different_modes/Omega_stack.csv");
V_cpp = load("../test_different_modes/V_stack.csv");
dot_Omega_cpp = load("../test_different_modes/dot_Omega_stack.csv");
dot_V_cpp = load("../test_different_modes/dot_V_stack.csv");

N_cpp = load("../test_different_modes/N_stack.csv");
C_cpp = load("../test_different_modes/C_stack.csv");
Qa_cpp = load("../test_different_modes/Qa_stack.csv");



error_Q = Q_cpp - QX;
error_Q_norm = norm(error_Q);
% display(Q_cpp)
% display(QX)
display(error_Q_norm)


error_r = r_cpp - rX;
error_r_norm = norm(error_r);
display(error_r_norm)


error_Omega = Omega_cpp - eta_X(1:3, :);
error_Omega_norm = norm(error_Omega);
display(error_Omega_norm)


error_V = V_cpp - eta_X(4:6, :);
error_V_norm = norm(error_V);
display(error_V_norm)



error_dot_Omega = dot_Omega_cpp - eta_dot_X(1:3, :);
error_dot_Omega_norm = norm(error_dot_Omega);
display(error_dot_Omega_norm)


error_dot_V = dot_V_cpp - eta_dot_X(4:6, :);
error_dot_V_norm = norm(error_dot_V);
display(error_dot_V_norm)




error_C = C_cpp - Lambda_X(1:3, :);
error_C_norm = norm(error_C);
display(error_C_norm)


error_N = N_cpp - Lambda_X(4:6, :);
error_N_norm = norm(error_N);
display(error_N_norm)


error_Qa = Qa_cpp - Qa_X;
error_Qa_norm = norm(error_Qa);
display(error_Qa_norm)






%%  Now check the variational part

error_Delta_rotations_norms_stack = [];
error_Delta_positions_norms_stack = [];
error_Delta_Omega_norms_stack     = [];
error_Delta_V_norms_stack         = [];
error_Delta_dot_Omega_norms_stack = [];
error_Delta_dot_V_norms_stack     = [];
error_Delta_N_norms_stack         = [];
error_Delta_C_norms_stack         = [];
error_Delta_Qa_norms_stack        = [];

for i=1:Const.dim_base
    
    Delta_rotations_cpp = load("../test_different_modes/Delta_rotation_stack_Delta" + (i-1) + ".csv");
    
    error_Delta_rotations = Delta_rotations_cpp - Delta_zeta_X(1:3, :, i);
    error_Delta_rotations_norm = norm(error_Delta_rotations);
    error_Delta_rotations_norms_stack = [error_Delta_rotations_norms_stack error_Delta_rotations_norm];




    Delta_positions_cpp = load("../test_different_modes/Delta_position_stack_Delta" + (i-1) + ".csv");

    error_Delta_positions = Delta_positions_cpp - Delta_zeta_X(4:6, :, i);
    error_Delta_positions_norm = norm(error_Delta_positions);
    error_Delta_positions_norms_stack = [error_Delta_positions_norms_stack error_Delta_positions_norm];






    Delta_Omega_cpp = load("../test_different_modes/Delta_Omega_stack_Delta" + (i-1) + ".csv");
    
    error_Delta_Omega = Delta_Omega_cpp - Delta_eta_X(1:3, :, i);
    error_Delta_Omega_norms_stack = [error_Delta_Omega_norms_stack norm(error_Delta_Omega)];




    Delta_V_cpp = load("../test_different_modes/Delta_V_stack_Delta" + (i-1) + ".csv");

    error_Delta_V = Delta_V_cpp - Delta_eta_X(4:6, :, i);
    error_Delta_V_norms_stack = [error_Delta_V_norms_stack norm(error_Delta_V)];







    Delta_dot_Omega_cpp = load("../test_different_modes/Delta_dot_Omega_stack_Delta" + (i-1) + ".csv");
    
    error_Delta_dot_Omega = Delta_dot_Omega_cpp - Delta_eta_dot_X(1:3, :, i);
    error_Delta_dot_Omega_norms_stack = [error_Delta_dot_Omega_norms_stack norm(error_Delta_dot_Omega)];




    Delta_dot_V_cpp = load("../test_different_modes/Delta_dot_V_stack_Delta" + (i-1) + ".csv");

    error_Delta_dot_V = Delta_dot_V_cpp - Delta_eta_dot_X(4:6, :, i);
    Delta_dot_V = Delta_eta_dot_X(4:6, :, i);
%     display(Delta_dot_V_cpp)
%     display(Delta_dot_V)
    error_Delta_dot_V_norms_stack = [error_Delta_dot_V_norms_stack norm(error_Delta_dot_V)];





    Delta_C_cpp = load("../test_different_modes/Delta_C_stack_Delta" + (i-1) + ".csv");
    
    error_Delta_C = Delta_C_cpp - Delta_Lambda_X(1:3, :, i);
    error_Delta_C_norms_stack = [error_Delta_C_norms_stack norm(error_Delta_C)];




    Delta_N_cpp = load("../test_different_modes/Delta_N_stack_Delta" + (i-1) + ".csv");

    error_Delta_N = Delta_N_cpp - Delta_Lambda_X(4:6, :, i);
    error_Delta_N_norms_stack = [error_Delta_N_norms_stack norm(error_Delta_N)];





    Delta_Qa_cpp = load("../test_different_modes/Delta_Qa_stack_Delta" + (i-1) + ".csv");

    error_Delta_Qa = Delta_Qa_cpp - Delta_Qa_X(:, :, i);
    error_Delta_Qa_norms_stack = [error_Delta_Qa_norms_stack norm(error_Delta_Qa)];

end

display(error_Delta_rotations_norms_stack)
display(error_Delta_positions_norms_stack)

display(error_Delta_Omega_norms_stack)
display(error_Delta_V_norms_stack)

display(error_Delta_dot_Omega_norms_stack)
display(error_Delta_dot_V_norms_stack)


display(error_Delta_C_norms_stack)
display(error_Delta_N_norms_stack)

display(error_Delta_Qa_norms_stack)





