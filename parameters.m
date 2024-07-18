
% Definition of parameters

phi_pl = 5;             % Payload fairing diameter [m]
L_pl = 12;              % Payload fairing length [m]
th_pl = 0.005;          % Payload fairing thickness [m]

M_pl = 7000;            % Payload mass [kg]

M_eng_vac = 500;        % Vacuum engine mass [kg]
M_eng_sl = 470;         % Sea level engine mass [kg]
L_eng_vac = 6;          % Longitude of vacuum engine [m]
L_eng_sl = 3;           % Longitude of sea level engine [m]
rho_prop = 1000;        % Density of the propellant [kg/m^3]

alpha = deg2rad(0.5);   % Angle of attack [deg->rad]
g0 = 9.81;              % Gravitational acceleration at Earth's surface [m/s^2]
thrust_1e = 900e3;      % Thrust of one SL engine [N]
SF = 2;                 % Safety factor

rho = 2550;             % Structural material density [kg/m^3]
E = 77e9;               % Young modulus [Pa]
sig_y = 210e6;          % Yield strength [Pa]
sig_s = 200e6;          % Shear strength [Pa]

h_tar = 300e3;          % Target orbital altitude [m]
v_tar = 7.73e3;         % Target orbital velocity [m/s]

