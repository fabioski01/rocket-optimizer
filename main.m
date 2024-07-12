clc, clear, close all

% Design varibles
phi = [3.7, 3.7, 3.7];    % Rocket stage diameter (3rd, 2nd, 1st) [m]
th = [0.005, 0.005, 0.005]; % Rocket stage wall thickness (3rd, 2nd, 1st) [m]
L = [10, 10, 20];     % Rocket stage length (3rd, 2nd, 1st) [m]
n_st = 3;   % Number of stages
n_e = 6;   % Number of engines (in the 1st stage)

% phi = [3.5, 3.5];    % Rocket stage diameter (3rd, 2nd, 1st) [m]
% th = [0.005, 0.005]; % Rocket stage wall thickness (3rd, 2nd, 1st) [m]
% L = [10, 30];     % Rocket stage length (3rd, 2nd, 1st) [m]
% n_st = 2;   % Number of stages
% n_e = 9;   % Number of engines (in the 1st stage)

% Define constant parameters
parameters;

% Objective function
f = structural_mass(phi,th,L,n_st,n_e);

% Constraints
[g,h] = constraints(phi,th,L,n_st,n_e);

% Sensitivity analysis
% gfd_mass(phi,th,L,n_st,n_e);

