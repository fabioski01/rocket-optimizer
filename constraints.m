function g = constraints(phi,th,L,n_st,n_e)

% Define constant parameters
parameters;

num = 100;              % Number of evaluated points per component
A = pi*(phi/2)^2;       % Reference area (rocket cross-section) [m^2]
thrust = thrust_1e*n_e; % Total SL thrust [N]

L_prime = L_pl*15/16;
l = [L_prime/3,   L_prime*2/3, L_prime/15,...
     repmat([L-L_eng_vac, L_eng_vac],1,n_st-1),...
     L-L_eng_sl,  L_eng_sl]; % Length of each component [m]

V1 = pi*phi_pl^2*l(1)/8 - pi*(phi_pl-2*th_pl)^2*(l(1)-th_pl)/8; % Volume of the nose cone of the fairing [m^3]
V2 = pi*(phi_pl/2)^2*l(2) - pi*(phi_pl/2 - th_pl)^2*l(2); % Volume of the cylindrical body of the fairing [m^3]
V3 = 1/3*pi*l(3)*((phi_pl/2)^2 + (phi/2)^2 + phi_pl/2*phi/2) - ...
     1/3*pi*l(3)*((phi_pl/2-th_pl)^2 + (phi/2-th_pl)^2 + (phi_pl/2-th_pl)*(phi/2-th_pl)); % Volume of the conical boattail of the fairing [m^3] 
V_tank = pi*(phi/2)^2*l(4) - pi*(phi/2 - th)^2*l(4); % Volume of the rocket shell for the lenght of a single tank [m^3]
V_eng_vac = pi*(phi/2)^2*l(5) - pi*(phi/2 - th)^2*l(5); % Volume of the rocket shell for the lenght of a vacuum engine [m^3]
V_eng_sl = pi*(phi/2)^2*l(end) - pi*(phi/2 - th)^2*l(end); % Volume of the rocket shell for the lenght of a sea level engine [m^3]

M = [V1*rho, V2*rho + M_pl, V3*rho,...
     repmat([V_tank*rho + M_prop/n_st, V_eng_vac*rho + M_eng_vac],1,n_st-1),...
     V_tank*rho + M_prop/n_st, V_eng_sl*rho + M_eng_sl*n_e]; % Mass of each component [kg]
m = M./l; % Mass of each component (per unit length) [kg/m]

cg_local = zeros(1,length(l)); % CoG of each component using local coord [m]
cg = zeros(1,length(l)); % CoG of each component using global coord [m]

Cp_local = zeros(1,length(l)); % Center of pressure of each component using local coord [m]
Cp = zeros(1,length(l)); % Center of pressure of each component using global coord [m]

Cna = zeros(1,length(l)); % Slope of the normal force coefficient at alpha = 0 of each component (per unit length) [1/m]
Cd = zeros(1,length(l)); % Drag coefficient of each component (per unit length) [1/m]

x1 = 0:0.1:l(1); theta = acos(1-2*x1/l(1)); 
y1 = phi_pl/2*sqrt(theta-sin(2*theta)/2)/sqrt(pi);
polyin1 = polyshape(x1,y1);
[cg_n1,~] = centroid(polyin1); % CoG of the nose cone of the fairing [m]

x3 = [0, 0, l(3), l(3)];
y3 = [0, phi_pl/2, phi/2, 0];
polyin3 = polyshape(x3,y3);
[cg_n3,~] = centroid(polyin3); % CoG of the conical boattail of the fairing [m]

suma = 0;
for i = 1:length(l)
    if (i == 1)
        cg_local(i) = cg_n1;
        cg(i) = cg_local(i) + suma;
        Cp_local(i) = l(i)/2;
        Cp(i) = Cp_local(i) + suma;
        Cna(i) = 2/l(i);
        Cd(i) = 0.1/l(i) + 0.15/sum(l);
    elseif (i == 3)
        cg_local(i) = cg_n3;
        cg(i) = cg_local(i) + suma;
        Cp_local(i) = l(i)/3*(1+(1-phi_pl/phi)/(1-(phi_pl/phi)^2));
        Cp(i) = Cp_local(i) + suma;
        Cna(i) = 8/(pi*phi_pl^2)*(pi*phi^2/4 - pi*phi_pl^2/4)/l(i);
        Cd(i) = 0.15/sum(l);
    else
        cg_local(i) = l(i)/2;
        cg(i) = cg_local(i) + suma;
        Cp_local(i) = l(i)/2;
        Cp(i) = Cp_local(i) + suma;
        Cna(i) = 0;
        Cd(i) = 0.15/sum(l);
    end
    suma = l(i) + suma;
end
CG = sum(m.*l.*cg)/sum(m.*l); % Total CoG [m]

t = zeros(1,length(l));
t(end) = thrust;
G = (t(end)-sum(q*A*Cd.*l)-sum(m.*l)*g0)/(sum(m.*l)*g0); % Axial acceleration [g's]

n = q*A*alpha*Cna;
acc_ang = (sum(n.*l)*sum(m.*l.*cg)/sum(m.*l)-sum(n.*l.*Cp))/(sum(m.*l.*(CG-cg))*sum(m.*l.*cg)/sum(m.*l)-sum(m.*l.*(CG-cg).*cg)); % Angular acceleration
Gy = (sum(n.*l)-acc_ang*sum(m.*l.*(CG-cg)))/(sum(m.*l)*g0); % Lateral acceleration [g's]


%% LOADS
x_last = 0;
D_last = 0;
Ix_last = 0;
F_last = 0;
Iy_last = 0;
N_last = 0;
for i = 1:length(l)
    x = linspace(0,l(i),num);
    pos(:,i) = (x + x_last).';

    D_tram = q*A*Cd(i)*x;
    Ix_tram = (G+1)*g0*m(i)*x;
    F_tram = t(i)*ones(1,num);
    Iy_tram = Gy*g0*m(i)*x + acc_ang*m(i)*x.*(CG-(x/2+x_last));
    col = find(x>=Cp_local(i));
    N_tram = n(i)*l(i)*[zeros(1,col(1)-1),ones(1,num-col(1)+1)];
    if (i>1)
        MIy_tram = Gy*g0*m(i)*x.*(x-x/2) + (sum(repmat(Gy*g0*m(1:i-1).*l(1:i-1),num,1).*((x+x_last).'-cg(1:i-1)),2)).'...
                 + acc_ang*m(i)*x.*(CG-(x/2+x_last)).*(x-x/2) + (sum(repmat(acc_ang*m(1:i-1).*l(1:i-1).*(CG-cg(1:i-1)),num,1).*((x+x_last).'-cg(1:i-1)),2)).';
        MN_tram = -N_tram.*(x-Cp_local(i)) - (sum(repmat(n(1:i-1).*l(1:i-1),num,1).*((x+x_last).'-Cp(1:i-1)),2)).';
    else
        MIy_tram = Gy*g0*m(i)*x.*(x-x/2) + acc_ang*m(i)*x.*(CG-(x/2+x_last)).*(x-x/2);
        MN_tram = -N_tram.*(x-Cp_local(i));
    end
    
    % Axial forces
    D(:,i) = (D_tram + D_last).';
    Ix(:,i) = (Ix_tram + Ix_last).';
    F(:,i) = (F_tram + F_last).';
    % Shear forces
    Iy(:,i) = (Iy_tram + Iy_last).';
    N(:,i) = (N_tram + N_last).';
    % Bending moments
    MIy(:,i) = (MIy_tram).';
    MN(:,i) = (MN_tram).';

    x_last = pos(end,i);
    D_last = D(end,i);
    Ix_last = Ix(end,i);
    F_last = F(end,i);
    Iy_last = Iy(end,i);
    N_last = N(end,i);
end

pos = reshape(pos,1,[]);
D = reshape(D,1,[]);
Ix = reshape(Ix,1,[]);
F = reshape(F,1,[]);
Iy = reshape(Iy,1,[]);
N = reshape(N,1,[]);
MIy = reshape(MIy,1,[]);
MN = reshape(MN,1,[]);

Axial = F - D - Ix;
figure
plot(pos,Axial/1e3)
hold on
plot(pos,zeros(1,length(pos)),'k-')
xlabel('x (m)','Interpreter','latex')
ylabel('Axial force (kN)','Interpreter','latex')
xlim([0 sum(l)])

Shear = N - Iy;
figure
plot(pos,Shear)
hold on
plot(pos,zeros(1,length(pos)),'k-')
xlabel('x (m)','Interpreter','latex')
ylabel('Shear force (N)','Interpreter','latex')
xlim([0 sum(l)])

Bending = MIy + MN;
figure
plot(pos,Bending/1e3)
hold on
plot(pos,zeros(1,length(pos)),'k--')
xlabel('x (m)','Interpreter','latex')
ylabel('Bending moment (kNm)','Interpreter','latex')
xlim([0 sum(l)])


%% STRESSES
Ac = pi*((phi/2)^2 - (phi/2-th)^2);
Z = pi*((phi/2)^4 - (phi/2-th)^4)/(4*phi/2);
sig_axial = Axial/Ac;
sig_shear = Shear/Ac;
sig_bending = Bending/Z;

I = pi/64*(phi^4 - (phi-2*th)^4); % Area moment of inertia of the cross-section [m^4]
sig_crit = pi*E*I/(L^2*phi*th); % Critical buckling stress

sig_c_max = abs(min(sig_axial+sig_bending));
sig_t_max = abs(max(sig_axial+sig_bending));
sig_s_max = abs(max(sig_shear));
fprintf('Maximum compressive stress (MPa): %.4f\n',sig_c_max/1e6)
fprintf('Maximum tensile stress (MPa): %.4f\n',sig_t_max/1e6)
fprintf('Maximum shear stress (MPa): %.4f\n',sig_s_max/1e6)


%% GEOMETRY
x2 = [x1(end), x1(end)+l(2)]; y2 = [y1(end), y1(end)];
x3 = [x2(end), x2(end)+l(3)]; y3 = [y1(end), phi/2];
x4 = [x3(end), x3(end)+sum(l(4:end)), x3(end)+sum(l(4:end))]; y4 = [y3(end), y3(end), 0];

figure
plot(x1,y1,'k')
hold on
plot(x2,y2,'k')
plot(x3,y3,'k')
plot(x4,y4,'k')
plot(x1,-y1,'k')
plot(x2,-y2,'k')
plot(x3,-y3,'k')
plot(x4,-y4,'k')
plot([x3(end), x3(end)],[-phi/2, phi/2],'k')
plot([sum(l(1:5)), sum(l(1:5))],[-phi/2, phi/2],'k')
plot([sum(l(1:7)), sum(l(1:7))],[-phi/2, phi/2],'k')
axis equal


%% CONSTRAINTS
g1 = 1 - sig_y/(SF*sig_c_max);
g2 = 1 - sig_y/(SF*sig_t_max);
g3 = 1 - sig_s/(SF*sig_s_max);
g4 = 1 - sig_crit/(SF*sig_c_max);
g5 = 20*th/phi - 1;
g6 = phi/L - 1;

g = [g1,g2,g3,g4,g5,g6];

end
