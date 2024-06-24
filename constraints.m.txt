clc, clear, close all

phi = 6; % rocket diameter
th = 0.005; % rocket wall thickness
l = [3,5,2,10,2];
m = [150,500,50,3000,80]./l;

% Variables
num = 100;
q = 0.5*1.225*100^2;
alpha = 0.5*pi/180;
A = pi*(phi/2)^2;
g0 = 9.81;

cg_local = [1.5,2.5,1,5,1];
cg = [1.5,2.5+3,1+3+5,5+3+5+2,1+3+5+2+10];
CG = sum(m.*l.*cg)/sum(m.*l);

Cd = [0.3,0.005,0.005,0.005,0.5];
Cna = [0.07,0.01,0.01,0.01,0.05];

Cp = [2,2.5+3,1+3+5,5+3+5+2,1+3+5+2+10];
Cp_local = [2,2.5,1,5,1];

t = [0,0,0,0,8e+05];
G = (t(end)-sum(q*A*Cd.*l)-sum(m.*l)*g0)/(sum(m.*l)*g0);

n = q*A*alpha*Cna;
acc_ang = (sum(n.*l)*sum(m.*l.*cg)/sum(m.*l)-sum(n.*l.*Cp))/(sum(m.*l.*(CG-cg))*sum(m.*l.*cg)/sum(m.*l)-sum(m.*l.*(CG-cg).*cg));
Gy = (sum(n.*l)-acc_ang*sum(m.*l.*(CG-cg)))/(sum(m.*l)*g0);

%% LOADS
x_last = 0;
D_last = 0;
I_last = 0;
F_last = 0;
Iy_last = 0;
N_last = 0;
for i = 1:length(l)
    x = linspace(0,l(i),num);
    pos(:,i) = (x + x_last).';

    D_tram = q*A*Cd(i)*x;
    I_tram = (G+1)*g0*m(i)*x;
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
    I(:,i) = (I_tram + I_last).';
    F(:,i) = (F_tram + F_last).';
    % Shear forces
    Iy(:,i) = (Iy_tram + Iy_last).';
    N(:,i) = (N_tram + N_last).';
    % Bending moments
    MIy(:,i) = (MIy_tram).';
    MN(:,i) = (MN_tram).';

    x_last = pos(end,i);
    D_last = D(end,i);
    I_last = I(end,i);
    F_last = F(end,i);
    Iy_last = Iy(end,i);
    N_last = N(end,i);
end

pos = reshape(pos,1,[]);
D = reshape(D,1,[]);
I = reshape(I,1,[]);
F = reshape(F,1,[]);
Iy = reshape(Iy,1,[]);
N = reshape(N,1,[]);
MIy = reshape(MIy,1,[]);
MN = reshape(MN,1,[]);

Axial = F - D - I;
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

sig_c_max = abs(min(sig_axial+sig_bending));
sig_t_max = abs(max(sig_axial+sig_bending));
sig_s_max = abs(max(sig_shear));
fprintf('Maximum compressive stress (MPa): %.4f\n',sig_c_max/1e6)
fprintf('Maximum tensile stress (MPa): %.4f\n',sig_t_max/1e6)
fprintf('Maximum shear stress (MPa): %.4f\n',sig_s_max/1e6)


