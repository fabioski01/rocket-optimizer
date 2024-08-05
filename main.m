clc, clear, close all

% Design variables
phi = [3.65, 3.7, 3.75];    % Rocket stage diameter (3rd, 2nd, 1st) [m]
th = [0.0045, 0.005, 0.0059]; % Rocket stage wall thickness (3rd, 2nd, 1st) [m]
L = [10, 15, 20];     % Rocket stage length (3rd, 2nd, 1st) [m]
n_st = 3;   % Number of stages
n_e = 9;   % Number of engines (in the 1st stage)
% kick_angle = deg2rad(9);    % Kick angle for starting the gravity turn [deg]
% kick_time = 30;     % Moment in time when the kick angle is induced to start gravity turn [s] (can also be a percentage of the time span)

% phi = [3.5, 3.5];    % Rocket stage diameter (3rd, 2nd, 1st) [m]
% th = [0.005, 0.005]; % Rocket stage wall thickness (3rd, 2nd, 1st) [m]
% L = [10, 30];     % Rocket stage length (3rd, 2nd, 1st) [m]
% n_st = 2;   % Number of stages
% n_e = 9;   % Number of engines (in the 1st stage)

% Define constant parameters
parameters;

% Objective function
% f = structural_mass(phi,th,L,n_st,n_e);

% Constraints
[g,h] = constraints(phi,th,L,n_st,n_e);

% Sensitivity analysis
% gfd_mass(phi,th,L,n_st,n_e);

% Function analysis
% range_phi = 12; % [m] range for variation of 1st stage diameter
% step_phi = 0.01*range_phi; % [m] step for variation of 1st stage diameter
% range_th = 0.05; % [m] range for variation of 1st stage thickness
% step_th = 0.01*range_th; % [m] step for variation of 1st stage thickness
% range_L = 50; % [m] range for variation of 1st stage length
% step_L = 0.01*range_L; % [m] step for variation of 1st stage length
% n_st_range = 3; % discrete number of stages (MAX IS 3)
% n_e_range = 15; % discrete number of engines 
% fun_analysis_mass(phi,range_phi, step_phi, th, range_th, step_th, L, range_L, step_L, n_st, n_st_range, n_e, n_e_range)
