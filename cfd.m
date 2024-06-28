function sens = cfd(pert,phi,th,L,n_st,n_e)
% cfd computes the central finite difference sensitivity for every
% perturbation given in 'pert' and returns it in 'sens'

% Design perturbations
h_phi = phi.'*pert;
h_th = th.'*pert;
h_L = L.'*pert;

% Compute cfd derivative
sens = zeros(3*n_st,length(pert));
for j = 1:n_st
    phi_pert1 = phi; phi_pert2 = phi;
    th_pert1 = th; th_pert2 = th;
    L_pert1 = L; L_pert2 = L;
    for i = 1:length(pert)
        phi_pert1(j) = phi(j) + h_phi(j,i);
        th_pert1(j) = th(j) + h_th(j,i);
        L_pert1(j) = L(j) + h_L(j,i);

        M_str11 = structural_mass(phi_pert1,th,L,n_st,n_e);
        M_str21 = structural_mass(phi,th_pert1,L,n_st,n_e);
        M_str31 = structural_mass(phi,th,L_pert1,n_st,n_e);

        phi_pert2(j) = phi(j) - h_phi(j,i);
        th_pert2(j) = th(j) - h_th(j,i);
        L_pert2(j) = L(j) - h_L(j,i);
    
        M_str12 = structural_mass(phi_pert2,th,L,n_st,n_e);
        M_str22 = structural_mass(phi,th_pert2,L,n_st,n_e);
        M_str32 = structural_mass(phi,th,L_pert2,n_st,n_e);
        
        sens(j,i) = (M_str11 - M_str12)/(2*h_phi(j,i));
        sens(j+n_st,i) = (M_str21 - M_str22)/(2*h_th(j,i));
        sens(j+n_st*2,i) = (M_str31 - M_str32)/(2*h_L(j,i));
    end
end

end