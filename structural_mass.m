function M_str = structural_mass(phi,th,L,n_st,n_e)

% Define parameters
parameters;

% Length of each fairing component [m]
L_prime = L_pl*15/16;
l = [L_prime/3, L_prime*2/3, L_prime/15];

% Volume of each fairing component [m^3]
% Nose cone
V1 = pi*phi_pl^2*l(1)/8 - pi*(phi_pl-2*th_pl)^2*(l(1)-th_pl)/8;
% Cylindrical body
V2 = pi*(phi_pl/2)^2*l(2) - pi*(phi_pl/2 - th_pl)^2*l(2);
% Conical boattail
V3 = 1/3*pi*l(3)*((phi_pl/2)^2 + (phi/2)^2 + phi_pl/2*phi/2) - ...
     1/3*pi*l(3)*((phi_pl/2-th_pl)^2 + (phi/2-th_pl)^2 + (phi_pl/2-th_pl)*(phi/2-th_pl));

% Mass of the fairing
M_f = rho*(V1+V2+V3);

% Compute structural mass
M_str = rho*L*pi*((phi/2)^2-(phi/2-th)^2)*n_st + M_f + M_eng_sl*n_e + M_eng_vac*(n_st-1);

end