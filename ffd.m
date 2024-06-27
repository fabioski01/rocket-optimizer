function sens = ffd(pert,phi,th,L,n_st,n_e)
% ffd computes the forward finite difference sensitivity for every
% perturbation given in 'pert' and returns it in 'sens'

% Design perturbations
h_phi = phi*pert;
h_th = th*pert;
h_L = L*pert;

% Compute structural mass for unperturbed case
M_str0 = structural_mass(phi,th,L,n_st,n_e);

% Compute ffd derivative
sens = zeros(3,length(pert));
for i = 1:length(pert)
    M_str1 = structural_mass(phi+h_phi(i),th,L,n_st,n_e);
    M_str2 = structural_mass(phi,th+h_th(i),L,n_st,n_e);
    M_str3 = structural_mass(phi,th,L+h_L(i),n_st,n_e);
    
    sens(1,i) = (M_str1 - M_str0)/(h_phi(i));
    sens(2,i) = (M_str2 - M_str0)/(h_th(i));
    sens(3,i) = (M_str3 - M_str0)/(h_L(i));
end

end