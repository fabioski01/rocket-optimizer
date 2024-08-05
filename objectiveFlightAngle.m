function cost = objectiveFlightAngle(x, phi,n_st,n_e,M_prop,M_str)
    % Extract optimized variables
    kick_time = x(1);
    kick_angle = x(2);

    % Definition of parameters
    
    phi_pl = 5;             % Payload fairing diameter [m]
    L_pl = 12;              % Payload fairing length [m]
    th_pl = 0.005;          % Payload fairing thickness [m]
    
    M_pl = 7000;            % Payload mass [kg]
    
    M_eng_vac = 500;        % Vacuum engine mass [kg]
    M_eng_sl = 470;         % Sea level engine mass [kg]
    L_eng_vac = 6;          % Longitude of vacuum engine [m]
    L_eng_sl = 3;           % Longitude of sea level engine [m]
    rho_prop = 1000;        % Density of the propellant [kg/m^3]
    
    alpha = deg2rad(0.5);   % Angle of attack [deg->rad]
    g0 = 9.81;              % Gravitational acceleration at Earth's surface [m/s^2]
    thrust_1e = 900e3;      % Thrust of one SL engine [N]
    SF = 2;                 % Safety factor
    
    rho = 2550;             % Structural material density [kg/m^3]
    E = 77e9;               % Young modulus [Pa]
    sig_y = 210e6;          % Yield strength [Pa]
    sig_s = 200e6;          % Shear strength [Pa]
    
    h_tar = 300e3;          % Target orbital altitude [m]
    v_tar = 7.73e3;         % Target orbital velocity [m/s]
    
    
    %% 1) Initial parameters and variables
    % Me
    Me = sym('Me','real');
    gamma = 1.1488;
    Ae_At = 21; % Exit to throat area ratio
    rel_arees = Ae_At == (2/(gamma+1))^(0.5*(gamma+1)/(gamma-1))*1/Me*(1+(gamma-1)/2*Me^2)^(0.5*(gamma+1)/(gamma-1));
    Me_V = vpasolve(rel_arees,Me);
    Me_V = double(Me_V); % Exit Mach number
    
    % CF & At
    pa_pc = 0.0094; % Initial pressure ratio
    CF_vac = (2/(gamma+1))^(0.5*(gamma+1)/(gamma-1))*(gamma*Me_V+1/Me_V)/sqrt(1+(gamma-1)/2*Me_V^2);
    MFP_Me = sqrt(gamma)*Me_V/(1+(gamma-1)/2*Me_V^2)^(0.5*(gamma+1)/(gamma-1));
    MFP_Mt = sqrt(gamma)*1/(1+(gamma-1)/2*1^2)^(0.5*(gamma+1)/(gamma-1));
    CF = CF_vac-pa_pc*MFP_Mt/MFP_Me;
    
    Wp = M_prop*g0; % Initial propellant weight [N]
    Ws = M_str*g0; % Initial structural weight [N]
    % Wu = 400*g0; % payload weight [N]
    Wu = M_pl*g0; % payload weight [N]
    W = sum(Wp)+sum(Ws)+Wu; % Initial weight [N]
    F = thrust_1e; % Initial thrust [N]
    Pc = 10.8e6; % Chamber pressure [Pa]
    At = F/(Pc*CF); % Throat area [m^2]
    Ae_sl = At*Ae_At;
    S = pi*(max(phi)/2)^2; % Rocket's maximum cross-section [m^2]
    
    % MASS FLOW RATE
    R = 8.314472; % Universal gas constant [J/molK]
    MM = 22.686; % Molecular mass [g/mol]
    Rg = R/(MM/1000); % Gas constant [J/kgK]
    Tc = 3655.73; % Chamber temperature [K]
    m_dot = MFP_Mt*Pc*At/sqrt(Rg*Tc);
    
    %% 2) Initial conditions
    h0 = 0.001; % [m]
    v0 = 0.001; % [m/s]
    ang0 = deg2rad(90); % [deg->rad]
    
    %% 3) Time span and time step
    t_span1 = Wp(end)/(m_dot*g0*n_e); % [s]
    t_step = 1e-3; % [s]

    %% 4) Solver
    % Solver - FIRST STAGE before kick
    t0 = 1e-6;
    [~, y_1st_initial] = ode45(@(t, y) Fsyst(t, t0, y, g0, W, m_dot, gamma, Rg, CF_vac, Pc, MFP_Me, MFP_Mt, At, S, n_e), 0:t_step:kick_time, [h0; v0; ang0]);
    
        % Extract final values after initial stage
    final_values_initial = y_1st_initial(end, :);  % Final [altitude, velocity, angle]
    h0_new = final_values_initial(1);
    v0_new = final_values_initial(2);
    ang0_new = ang0 - kick_angle;

    % Solver - FIRST STAGE after gravity turn
    t0 = kick_time;
    [~, y_1st_final] = ode45(@(t, y) Fsyst(t, t0, y, g0, W, m_dot, gamma, Rg, CF_vac, Pc, MFP_Me, MFP_Mt, At, S, n_e), kick_time:t_step:t_span1, [h0_new; v0_new; ang0_new]);

    % Solver - SECOND STAGE
        % Me
    Me = sym('Me','real');
    gamma = 1.1488;
    Ae_At = 116; % Exit to throat area ratio
    rel_arees = Ae_At == (2/(gamma+1))^(0.5*(gamma+1)/(gamma-1))*1/Me*(1+(gamma-1)/2*Me^2)^(0.5*(gamma+1)/(gamma-1));
    Me_V = vpasolve(rel_arees,Me);
    Me_V = double(Me_V); % Exit Mach number
    
        % CF & At
    pa_pc = 0.0094; % Initial pressure ratio
    CF_vac = (2/(gamma+1))^(0.5*(gamma+1)/(gamma-1))*(gamma*Me_V+1/Me_V)/sqrt(1+(gamma-1)/2*Me_V^2);
    MFP_Me = sqrt(gamma)*Me_V/(1+(gamma-1)/2*Me_V^2)^(0.5*(gamma+1)/(gamma-1));
    MFP_Mt = sqrt(gamma)*1/(1+(gamma-1)/2*1^2)^(0.5*(gamma+1)/(gamma-1));
    CF = CF_vac-pa_pc*MFP_Mt/MFP_Me;
    
    F = thrust_1e; % Initial thrust [N]
    Pc = 10.8e6; % Chamber pressure [Pa]
    At = F/(Pc*CF); % Throat area [m^2]
    Ae_vac = At*Ae_At;
    S = pi*(max(phi(1:end-1))/2)^2; % Rocket's maximum cross-section [m^2]
    
        % MASS FLOW RATE
    R = 8.314472; % Universal gas constant [J/molK]
    MM = 22.686; % Molecular mass [g/mol]
    Rg = R/(MM/1000); % Gas constant [J/kgK]
    Tc = 3655.73; % Chamber temperature [K]
    m_dot = MFP_Mt*Pc*At/sqrt(Rg*Tc);
    t_span23 = Wp(1:end-1)./(m_dot*g0); % [s]
    
    second_stage_initial = y_1st_final(end,:);  % Final [altitude, velocity, angle]
    % New initial conditions for the integration of second stage
    h0_2nd = second_stage_initial(1);
    v0_2nd = second_stage_initial(2);
    ang0_2nd = second_stage_initial(3);
    % Solve system of 2nd stage
    W = W - Ws(end) - Wp(end);
    t0 = t_span1;
    [~, y_2nd] = ode45(@(t,y) Fsyst(t,t0,y,g0,W,m_dot,gamma,Rg,CF_vac,Pc,MFP_Me,MFP_Mt,At,S,1), t_span1:t_step:(t_span1+t_span23(1)), [h0_2nd;v0_2nd;ang0_2nd]);
    
    if (n_st == 3)
        % THIRD STAGE
        third_stage_initial = y_2nd(end,:);  % Final [altitude, velocity, angle]
        % New initial conditions for the integration of second stage
        h0_3rd = third_stage_initial(1);
        v0_3rd = third_stage_initial(2);
        ang0_3rd = third_stage_initial(3);
        % Solve system of 3rd stage
        W = W - Ws(end-1) - Wp(end-1);
        t0 = (t_span1+t_span23(end-1));
        
        % Adjust ODE solver options % This ODE causes problems
        opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
        [~, y_3rd] = ode45(@(t,y) Fsyst(t,t0,y,g0,W,m_dot,gamma,Rg,CF_vac,Pc,MFP_Me,MFP_Mt,At,S,1), (t_span1+t_span23(1)):t_step:(t_span1+sum(t_span23)), [h0_3rd;v0_3rd;ang0_3rd], opts);
    
        inertial_stage_initial = y_3rd(end,:);
        h0_inr = inertial_stage_initial(1);
        v0_inr = inertial_stage_initial(2);
        ang0_inr = inertial_stage_initial(3);
        [~, y_inr] = ode45(@(t,y) Fsyst2(y), (t_span1+sum(t_span23)):t_step:(t_span1+sum(t_span23))+20, [h0_inr;v0_inr;ang0_inr]);
    end

    % Extract final values after first stage
    if (n_st == 3)
        final_values = y_3rd(end, :); % Final [altitude, velocity, angle]
    else
        final_values = y_2nd(end, :); % Final [altitude, velocity, angle]
    end

    % The cost is the absolute difference between the final flight path angle and zero
    final_flight_path_angle = final_values(3);
    cost = abs(final_flight_path_angle);
end