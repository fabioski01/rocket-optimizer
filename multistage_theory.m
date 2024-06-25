%%% Multi-Stage rocket theory
clc, clear, close all

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
Mp = 60000; % propellant mass in kg, only fuel and oxidizer

% constraints derived from 1st order ones, which are calculated depending on variables
Mass_ratio = exp(DeltaV/c_eff) % from ideal rocket equation M0/Me
phi = 1 - (1/Mass_ratio) % propellant ratio Mp/M0

% computations
M0 = Mc + Mp + Mu % gross lift off weight, to be minimized
Me = M0 - Mp % final mass
lambda = Mu/M0 % payload mass ratio >1, ideally high;
epsilon = Mc/(Mc+Mp) % construction mass ratio <1,  ideally low;
phi_check = Mp/M0 % propellant mass ratio <1. Phi = (1-epsilon)*(1-lambda)
Mass_ratio_check = M0/Me % total mass ratio



