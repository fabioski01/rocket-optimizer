clear; close all; clc;

%% AN ANALITICAL APPROACH OF A ROCKET TRAJECTORY
% Jan Canet, Marc Massons, Eduardo Jiménez, Albert Llonch

% 1) Initial parameters and variables
%% Me
Me = sym('Me','real');
gamma = 1.290475;
Ae_At = 30; % Exit to throat area ratio
rel_arees = Ae_At == (2/(gamma+1))^(0.5*(gamma+1)/(gamma-1))*1/Me*(1+(gamma-1)/2*Me^2)^(0.5*(gamma+1)/(gamma-1));
Me_V = vpasolve(rel_arees,Me);
Me_V = double(Me_V); % Exit Mach number

%% CF & At
pa_pc = 1/71631; % Initial pressure ratio
CF_vac = (2/(gamma+1))^(0.5*(gamma+1)/(gamma-1))*(gamma*Me_V+1/Me_V)/sqrt(1+(gamma-1)/2*Me_V^2);
MFP_Me = sqrt(gamma)*Me_V/(1+(gamma-1)/2*Me_V^2)^(0.5*(gamma+1)/(gamma-1));
MFP_Mt = sqrt(gamma)*1/(1+(gamma-1)/2*1^2)^(0.5*(gamma+1)/(gamma-1));
CF = CF_vac-pa_pc*MFP_Mt/MFP_Me;

g50 = gravity(50e3);
Wp = 50940*g50; % Initial propellant weight [N]
Ws = 37017*g50; % Initial structural weight [N]
W = Wp+Ws; % Initial weight [N]
F = W*1.25; % Initial thrust [N]
Pc = 5.44e6; % Chamber pressure [Pa]
At = F/(Pc*CF); % Throat area [m^2]
S = pi*(2.95/2)^2; % Rocket's maximum cross-section [m^2]

%% MASS FLOW RATE
R = 8.314472; % Universal gas constant [J/molK]
MM = 22.774; % Molecular mass [g/mol]
Rg = R/(MM/1000); % Gas constant [J/kgK]
Tc = 3594.6; % Chamber temperature [K]
m_dot = MFP_Mt*Pc*At/sqrt(Rg*Tc);

% 2) Initial conditions

h0 = 50e3; % [m]
v0 = h0/118; % [m/s]
ang0 = deg2rad(89); % [deg->rad]

% 3) Time span and time step

t_span = Wp/(m_dot*g50); % [s]
t_step = 1e-3; % [s]

% 4) Solver

[t0,y0] = ode45(@(t,y) Fsyst(t,y,g50,W,m_dot,gamma,Rg,CF_vac,Pc,MFP_Me,MFP_Mt,At,S),0:t_step:t_span,[h0;v0;ang0]);

% 5) Postprocess

figure
subplot(1,3,1)
plot(t0,y0(:,1)/1000)
xlim([0 t_span])
title('Altitude over time')
xlabel('Time [s]')
ylabel('Altitude [km]')
grid on

subplot(1,3,2)
plot(t0,y0(:,2))
xlim([0 t_span])
title('Velocity over time')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on

subplot(1,3,3)
plot(t0,y0(:,3)*180/pi)
xlim([0 t_span])
title('Flight path angle over time')
xlabel('Time [s]')
ylabel('Flight path angle [º]')
grid on

save t0.mat
save y0.mat










