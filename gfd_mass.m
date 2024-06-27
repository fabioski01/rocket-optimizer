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
Effd_phi = 100*(sens_forward(1,:)-dM_dphi)/dM_dphi;
Effd_th = 100*(sens_forward(2,:)-dM_dth)/dM_dth;
Effd_L = 100*(sens_forward(3,:)-dM_dL)/dM_dL;

Ecfd_phi = 100*(sens_central(1,:)-dM_dphi)/dM_dphi;
Ecfd_th = 100*(sens_central(2,:)-dM_dth)/dM_dth;
Ecfd_L = 100*(sens_central(3,:)-dM_dL)/dM_dL;

% Plot the results
% Sensitivity wrt phi
figure('color','w')
semilogx(pert, Effd_phi,'ro--');
hold on
semilogx(pert, Ecfd_phi,'b*--');
% Some additional lines
semilogx(pert, zeros(size(pert)),'k--');
semilogx(pert, ones(size(pert)),'k--');
semilogx(pert, -ones(size(pert)),'k--');

title(['Relative error of GFD sensitivities for mesh with ',...
    num2str(n_st),' stages and ',num2str(n_e),' engines.']);
xlabel('Relative design perturbation');
ylabel('Error [%] in dM/dphi')
legend('Forward FD','Central FD');

set(gca,'ylim',[-5 5]);

% Sensitivity wrt th
figure('color','w')
semilogx(pert, Effd_th,'ro--');
hold on
semilogx(pert, Ecfd_th,'b*--');
% Some additional lines
semilogx(pert, zeros(size(pert)),'k--');
semilogx(pert, ones(size(pert)),'k--');
semilogx(pert, -ones(size(pert)),'k--');

title(['Relative error of GFD sensitivities for mesh with ',...
    num2str(n_st),' stages and ',num2str(n_e),' engines.']);
xlabel('Relative design perturbation');
ylabel('Error [%] in dM/dth')
legend('Forward FD','Central FD');

set(gca,'ylim',[-5 5]);

% Sensitivity wrt L
figure('color','w')
semilogx(pert, Effd_L,'ro--');
hold on
semilogx(pert, Ecfd_L,'b*--');
% Some additional lines
semilogx(pert, zeros(size(pert)),'k--');
semilogx(pert, ones(size(pert)),'k--');
semilogx(pert, -ones(size(pert)),'k--');

title(['Relative error of GFD sensitivities for mesh with ',...
    num2str(n_st),' stages and ',num2str(n_e),' engines.']);
xlabel('Relative design perturbation');
ylabel('Error [%] in dM/dL')
legend('Forward FD','Central FD');

set(gca,'ylim',[-5 5]);

end