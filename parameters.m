
% Definition of parameters

phi_pl = 5;             % Payload fairing diameter [m]
L_pl = 12;              % Payload fairing length [m]
th_pl = 0.005;          % Payload fairing thickness [m]

M_pl = 15000;           % Payload mass [kg]
% M_prop = 300000;        % Propellant mass [kg]
M_prop = [100000, 100000, 100000];        % Propellant mass [kg]
M_eng_vac = 600;        % Vacuum engine mass [kg]
M_eng_sl = 150;         % Sea level engine mass [kg]

L_eng_vac = 3;          % Longitude of vacuum engine [m]
L_eng_sl = 1.5;         % Longitude of sea level engine [m]

q = 0.5*1.225*100^2;    % Maximum dynamic pressure [Pa]
alpha = deg2rad(0.5);   % Angle of attack [deg->rad]
g0 = 9.81;              % Gravitational acceleration at Earth's surface [m/s^2]
thrust_1e = 90e3;       % Thrust of one SL engine [N]
SF = 2;                 % Safety factor

rho = 2550;             % Structural material density [kg/m^3]
E = 77e9;               % Young modulus [Pa]
sig_y = 210e6;          % Yield strength [Pa]
sig_s = 200e6;          % Shear strength [Pa]

