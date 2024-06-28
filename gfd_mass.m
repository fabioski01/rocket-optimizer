function gfd_mass(phi,th,L,n_st,n_e)

% Use of logarithmic sensitivities?
% Complex perturbations needed?

% Define constant parameters
parameters;

% Analytical sensitivity values (used as reference)
dM_dphi = rho*L.*pi.*th;
dM_dth = rho*L.*pi.*(phi-2*th);
dM_dL = rho*pi*((phi/2).^2 -(phi/2-th).^2);

% Relative design perturbations
pert = 10.^[-15:-1];

% Compute forward finite difference sensitivities
disp('Forward finite differences')
tic
sens_forward = ffd(pert,phi,th,L,n_st,n_e);
toc

% Compute central finite difference sensitivities
disp('Central finite differences')
tic
sens_central = cfd(pert,phi,th,L,n_st,n_e);
toc

% Compute relative error [%]
for j = 1:n_st
    Effd_phi(j,:) = 100*(sens_forward(j,:)-dM_dphi(j))/dM_dphi(j);
    Effd_th(j,:) = 100*(sens_forward(j+n_st,:)-dM_dth(j))/dM_dth(j);
    Effd_L(j,:) = 100*(sens_forward(j+n_st*2,:)-dM_dL(j))/dM_dL(j);
    
    Ecfd_phi(j,:) = 100*(sens_central(j,:)-dM_dphi(j))/dM_dphi(j);
    Ecfd_th(j,:) = 100*(sens_central(j+n_st,:)-dM_dth(j))/dM_dth(j);
    Ecfd_L(j,:) = 100*(sens_central(j+n_st*2,:)-dM_dL(j))/dM_dL(j);
end

% Plot the results
for j = 1:n_st
    % Sensitivity wrt phi
    figure('color','w')
    semilogx(pert, Effd_phi(j,:),'ro--');
    hold on
    semilogx(pert, Ecfd_phi(j,:),'b*--');
    % Some additional lines
    semilogx(pert, zeros(size(pert)),'k--');
    semilogx(pert, ones(size(pert)),'k--');
    semilogx(pert, -ones(size(pert)),'k--');
    
    title(['Relative error of GFD sensitivities for a rocket with ',...
        num2str(n_st),' stages and ',num2str(n_e),' engines.']);
    subtitle(['Stage ',num2str(n_st-(j-1))])
    xlabel('Relative design perturbation');
    ylabel('Error [%] in dM/dphi')
    legend('Forward FD','Central FD');
    
    set(gca,'ylim',[-5 5]);
    
    % Sensitivity wrt th
    figure('color','w')
    semilogx(pert, Effd_th(j,:),'ro--');
    hold on
    semilogx(pert, Ecfd_th(j,:),'b*--');
    % Some additional lines
    semilogx(pert, zeros(size(pert)),'k--');
    semilogx(pert, ones(size(pert)),'k--');
    semilogx(pert, -ones(size(pert)),'k--');
    
    title(['Relative error of GFD sensitivities for a rocket with ',...
        num2str(n_st),' stages and ',num2str(n_e),' engines.']);
    subtitle(['Stage ',num2str(n_st-(j-1))])
    xlabel('Relative design perturbation');
    ylabel('Error [%] in dM/dth')
    legend('Forward FD','Central FD');
    
    set(gca,'ylim',[-5 5]);
    
    % Sensitivity wrt L
    figure('color','w')
    semilogx(pert, Effd_L(j,:),'ro--');
    hold on
    semilogx(pert, Ecfd_L(j,:),'b*--');
    % Some additional lines
    semilogx(pert, zeros(size(pert)),'k--');
    semilogx(pert, ones(size(pert)),'k--');
    semilogx(pert, -ones(size(pert)),'k--');
    
    title(['Relative error of GFD sensitivities for a rocket with ',...
        num2str(n_st),' stages and ',num2str(n_e),' engines.']);
    subtitle(['Stage ',num2str(n_st-(j-1))])
    xlabel('Relative design perturbation');
    ylabel('Error [%] in dM/dL')
    legend('Forward FD','Central FD');
    
    set(gca,'ylim',[-5 5]);
end

end