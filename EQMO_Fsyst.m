function Z = Fsyst(t,y,g50,W,m_dot,gamma,Rg,CF_vac,Pc,MFP_Me,MFP_Mt,At,S)
% Function that evaluates the two differential equations

Z = zeros(2,1);

m0 = W/g50;
m = m0-m_dot*t;
[Pext,rho_ext,Text] = Compute_P_rho_Text(y(1));
Pa = Pext;
M = y(2)/sqrt(gamma*Rg*Text);
Cd = ComputeCd(M);
g = gravity(y(1));

Z(1) = y(2);

Z(2) = 1/m*((CF_vac-Pa/Pc*MFP_Mt/MFP_Me)*Pc*At)-rho_ext*y(2)^2*S*Cd/(2*m)-g;

end

