function gfd_mass(n_st,n_e)

% Use of logarithmic sensitivities?
% Complex perturbations needed?

% Continuous design varibles
phi = 4;    % Rocket stage diameter [m]
th = 0.005; % Rocket stage wall thickness [m]
L = 20;     % Rocket stage length [m]

% Define constant parameters
parameters;

% Analytical sensitivity values (used as reference)
dM_dphi = rho*L*pi*th*n_st;
dM_dth = rho*L*pi*(phi-2*th)*n_st;
dM_dL = rho*pi*((phi/2)^2 -(phi/2-th)^2)*n_st;

% Relative design perturbations
pert = 10.^[-15:-1];

% Design perturbations
h_phi = phi*pert;
h_th = th*pert;
h_L = L*pert;

% Compute forward finite difference sensitivities:
disp('Forward finite differences')
tic
sens_forward = ffd();
toc

% Compute central finite difference sensitivities:
disp('Central finite differences')
tic
sens_central = cfd();
toc

% M_str = structural_mass(phi,th,L,n_st,n_e);

end