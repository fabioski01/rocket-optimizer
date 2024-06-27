%%% Multi-Stage rocket theory
clc, clear, close all
format long
% Goal: Minimization of gross lift-off weight

% 1st order Constraint: Specified mission constraints
Mu = 400; % Useful Payload mass in kg
h_t = 300; % Target orbit altitude, assuming circular orbit, in km
Ve = 0; % ideal end velocity (gravity-free space, vacuum) in km/s

% computations from constraints
GM = 3.986004418E14; % gravitational parameter Earth in m^3/s^2
Re = 6371.0; % Earth mean radius in km
r_orbit = (Re + h_t)*1000; % circular orbit radius in m
Vc = sqrt(GM/r_orbit)/1000; % needed orbit circular velocity in km/s
DeltaV = Vc - Ve; % needed deltaV for specified orbit

% initializing variables
c_eff = 4; % effective velocity in km/s depends on chosen propulsion system, example for LH2-LOX 4km/s is typical
Mc = 10000; % construction mass in kg, includes structural and engines but not propellant
% Mp = 60000; % propellant mass in kg, only fuel and oxidizer

% constraints derived from 1st order ones, which are calculated depending on variables
Mass_ratio = exp(DeltaV/c_eff) % from ideal rocket equation M0/Me
phi = 1 - (1/Mass_ratio) % propellant ratio Mp/M0
 
% % computations
% M0 = Mc + Mp + Mu % gross lift off weight, to be minimized
% Me = M0 - Mp % final mass
% lambda = Mu/M0 % payload mass ratio >1, ideally high;
% epsilon = Mc/(Mc+Mp) % construction mass ratio <1,  ideally low;
% phi_check = Mp/M0 % propellant mass ratio <1. Phi = (1-epsilon)*(1-lambda)
% Mass_ratio_check = M0/Me % total mass ratio

% Define the objective function
objective_function = @(Mp) calculate_mass_ratio_difference(Mp, Mc, Mu, Mass_ratio);

% Set initial guess and bounds for Mp
Mp_initial = 60000; % initial guess for propellant mass in kg
Mp_lower_bound = 0; % lower bound for propellant mass
Mp_upper_bound = 1e6; % upper bound for propellant mass

% Use fmincon to find the optimal Mp
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'TolFun', 1e-9, 'TolX', 1e-9, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 3000);
Mp_optimized = fmincon(objective_function, Mp_initial, [], [], [], [], Mp_lower_bound, Mp_upper_bound, [], options);

% Display the optimized Mp
fprintf('Optimized Mp: %.2f kg\n', Mp_optimized);

% Recalculate final values with optimized Mp
[M0, Me, lambda, ~, phi_check, Mass_ratio_check] = calculate_final_values(Mp_optimized, Mc, Mu);

fprintf('Mass_ratio: %.6f\n', Mass_ratio);
fprintf('Mass_ratio_check: %.6f\n', Mass_ratio_check);

% Assuming N=3 (3 stages)
% Assuming identical stages (same c_eff, epsilon, lambda, phi, Mass_ratio)
N=3; % number of stages
[M0_tot, Me_tot, lambda_tot, ~, ~, ~] = calculate_final_values(Mp_optimized, Mc, Mu) % for entire rocket
lambda_stage = lambda_tot^(1/N) % because lambda_tot equal the product of the lambda of each stage (for identical stages)

% other calculation for stage ratios
DeltaV_stage = DeltaV / N;
Mass_ratio_stage = exp(DeltaV_stage/c_eff)
phi_stage =  1-(1/Mass_ratio_stage)
epsilon_stage = 1 - phi_stage/(1-lambda_stage)
% lambda_stage_check = (exp(-Vc/(N*c_eff))-epsilon_stage) / (1-epsilon_stage)

% third stage
M0_third = Mu/lambda_stage % initial mass of the last (third) stage
% phi_stage = phi; % identical stages
% Mp_third = phi_stage*M0_third % propellant mass third stage
% epsilon_stage = epsilon_tot; % identical stages
% Mc_third = Mp_third*epsilon_stage % construction mass third stage
[Mp_third, Mc_third] = calculate_stage_values(phi_stage, M0_third, epsilon_stage)
% M0_third_check = Mc_third + Mp_third + Mu % check of M0 stage 

% second stage
M0_second = M0_third/lambda_stage % initial mass of the second stage
[Mp_second, Mc_second] = calculate_stage_values(phi_stage, M0_second, epsilon_stage)
% M0_second_check = Mc_second + Mp_second + M0_third % check of M0 stage 

% first stage
M0_first = M0_second/lambda_stage % initial mass of the second stage
[Mp_first, Mc_first] = calculate_stage_values(phi_stage, M0_first, epsilon_stage)
% M0_first_check = Mc_first + Mp_first + M0_second % check of M0 stage 

% Objective function to calculate the difference between Mass_ratio_check and Mass_ratio
function diff = calculate_mass_ratio_difference(Mp, Mc, Mu, Mass_ratio)
    [M0, Me, ~, ~, ~, Mass_ratio_check] = calculate_final_values(Mp, Mc, Mu);
    diff = abs(Mass_ratio_check - Mass_ratio);
end

% Function to calculate final values
function [M0, Me, lambda, epsilon, phi_check, Mass_ratio_check] = calculate_final_values(Mp, Mc, Mu)
    M0 = Mc + Mp + Mu; % gross lift off weight, to be minimized
    Me = M0 - Mp; % final mass
    lambda = Mu/M0; % payload mass ratio <1, ideally high;
    epsilon = Mc/(Mc+Mp); % construction mass ratio <1, ideally low;
    phi_check = Mp/M0; % propellant mass ratio <1. Phi = (1-epsilon)*(1-lambda)
    Mass_ratio_check = M0/Me; % total mass ratio
end

% Function for stage values
function [Mp_stage, Mc_stage] = calculate_stage_values(phi, M0_stage, epsilon)
    phi_stage = phi; % identical stages
    Mp_stage = phi_stage*M0_stage; % propellant mass stage
    epsilon_stage = epsilon; % identical stages
    Mc_stage = Mp_stage*epsilon_stage /( 1-epsilon_stage); % construction mass stage
end

