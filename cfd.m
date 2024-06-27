function sens = cfd(pert,phi,th,L,n_st,n_e)
% cfd computes the central finite difference sensitivity for every
% perturbation given in 'pert' and returns it in 'sens'

% Design perturbations
h_phi = phi*pert;
h_th = th*pert;
h_L = L*pert;

% Compute cfd derivative
sens = zeros(3,length(pert));
for i = 1:length(pert)
    M_str11 = structural_mass(phi+h_phi(i),th,L,n_st,n_e);
    M_str21 = structural_mass(phi,th+h_th(i),L,n_st,n_e);
    M_str31 = structural_mass(phi,th,L+h_L(i),n_st,n_e);

    M_str12 = structural_mass(phi-h_phi(i),th,L,n_st,n_e);
    M_str22 = structural_mass(phi,th-h_th(i),L,n_st,n_e);
    M_str32 = structural_mass(phi,th,L-h_L(i),n_st,n_e);
    
    sens(1,i) = (M_str11 - M_str12)/(2*h_phi(i));
    sens(2,i) = (M_str21 - M_str22)/(2*h_th(i));
    sens(3,i) = (M_str31 - M_str32)/(2*h_L(i));
end

end