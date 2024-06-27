clc, clear, close all

% Design varibles
% phi = 4;    % Rocket stage diameter [m]
% th = 0.005; % Rocket stage wall thickness [m]
% L = 20;     % Rocket stage length [m]

phi = [4, 5, 6];    % Rocket stage diameter [m]
th = [0.005, 0.005, 0.005]; % Rocket stage wall thickness [m]
L = [20, 20, 20];     % Rocket stage length [m]
n_st = 3;   % Number of stages
n_e = 9;    % Number of engines (in the 1st stage)

% Define constant parameters
parameters;

% Objective function
% f = structural_mass(phi,th,L,n_st,n_e);

% Constraints
g = constraints(phi,th,L,n_st,n_e);



