function Z = Fsyst_gravity(t,y,g50,W,m_dot,gamma,Rg,CF_vac,Pc,MFP_Me,MFP_Mt,At,S)
% Function that evaluates the two differential equations

% Z = zeros(2,1);
% 
% m0 = W/g50;
% m = m0-m_dot*t;
% [Pext,rho_ext,Text] = Compute_P_rho_Text(y(1));
% Pa = Pext;
% M = y(2)/sqrt(gamma*Rg*Text);
% Cd = ComputeCd(M);
% g = gravity(y(1));
% 
% Z(1) = y(2);
% 
% Z(2) = 1/m*((CF_vac-Pa/Pc*MFP_Mt/MFP_Me)*Pc*At)-rho_ext*y(2)^2*S*Cd/(2*m)-g;

Z = zeros(3,1);

Rt = 6371e3; % [m]

m0 = W/g50;
m = m0-m_dot*t;
[Pext,rho_ext,Text] = Compute_P_rho_Text(y(1));
Pa = Pext;
M = y(2)/sqrt(gamma*Rg*Text);
Cd = ComputeCd(M);
g = gravity(y(1));

% Altitude rate (vertical velocity component)
Z(1) = y(2)*sin(y(3));

% Velocity rate
% Thrust = (CF_vac - Pa / Pc * MFP_Mt / MFP_Me) * Pc * At;
% Drag = 0.5 * rho_ext * y(2)^2 * S * Cd;
% Z(2) = Thrust / m - Drag / m - g * sin(y(3));
Z(2) = 1/m*((CF_vac-Pa/Pc*MFP_Mt/MFP_Me)*Pc*At)-rho_ext*y(2)^2*S*Cd/(2*m)-g*sin(y(3));

% Flight path angle rate (gravity turn)
Z(3) = -(g/y(2) - y(2)/(Rt+y(1)))*cos(y(3));

end

