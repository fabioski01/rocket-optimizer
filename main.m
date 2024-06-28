clc, clear, close all

% Design varibles
phi = [4, 5, 6];    % Rocket stage diameter (3rd, 2nd, 1st) [m]
th = [0.005, 0.005, 0.005]; % Rocket stage wall thickness (3rd, 2nd, 1st) [m]
L = [20, 20, 20];     % Rocket stage length (3rd, 2nd, 1st) [m]
n_st = 3;   % Number of stages
n_e = 9;    % Number of engines (in the 1st stage)

% Define constant parameters
parameters;

% Objective function
% f = structural_mass(phi,th,L,n_st,n_e);

% Constraints
% g = constraints(phi,th,L,n_st,n_e);

% Sensitivity analysis
gfd_mass(phi,th,L,n_st,n_e);

