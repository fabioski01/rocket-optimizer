function sens = ffd(pert,phi,th,L,n_st,n_e)
% ffd computes the forward finite difference sensitivity for every
% perturbation given in 'pert' and returns it in 'sens'

% Design perturbations
h_phi = phi.'*pert;
h_th = th.'*pert;
h_L = L.'*pert;

% Compute structural mass for unperturbed case
M_str0 = structural_mass(phi,th,L,n_st,n_e);

% Compute ffd derivative
sens = zeros(3*n_st,length(pert));
for j = 1:n_st
    phi_pert = phi;
    th_pert = th;
    L_pert = L;
    for i = 1:length(pert)
        phi_pert(j) = phi(j) + h_phi(j,i);
        th_pert(j) = th(j) + h_th(j,i);
        L_pert(j) = L(j) + h_L(j,i);

        M_str1 = structural_mass(phi_pert,th,L,n_st,n_e);
        M_str2 = structural_mass(phi,th_pert,L,n_st,n_e);
        M_str3 = structural_mass(phi,th,L_pert,n_st,n_e);
        
        sens(j,i) = (M_str1 - M_str0)/(h_phi(j,i));
        sens(j+n_st,i) = (M_str2 - M_str0)/(h_th(j,i));
        sens(j+n_st*2,i) = (M_str3 - M_str0)/(h_L(j,i));
    end
end

end