function Z = Fsyst(t,t0,y,g0,W,m_dot,gamma,Rg,CF_vac,Pc,MFP_Me,MFP_Mt,At,S,n_e)
% Function that evaluates the two differential equations

Z = zeros(3,1);

Rt = 6371e3; % [m]

m0 = W/g0;
m = m0-m_dot*n_e*(t-t0);
[Pext,rho_ext,Text] = Compute_P_rho_Text(y(1));
Pa = Pext;
M = y(2)/sqrt(gamma*Rg*Text);
Cd = ComputeCd(M);
g = gravity(y(1));

Z(1) = y(2)*sin(y(3));

Z(2) = 1/m*((CF_vac-Pa/Pc*MFP_Mt/MFP_Me)*Pc*At*n_e)-rho_ext*y(2)^2*S*Cd/(2*m)-g*sin(y(3));

Z(3) = -(g/y(2) - y(2)/(Rt+y(1)))*cos(y(3));

end

