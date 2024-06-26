function [Pext, rho, Text] = Compute_P_rho_Text(z)
Rg = 8.31;
MM = 28.9644/1000; % [kg/mol]
r = Rg/MM;
z = z/1000; % eqs works in km

% if z>=0 && z<11
if z<11 % to keep going even with negative altitude errors
    Pext = 101325*(288.15/(288.15-6.5*z))^(34.1632/(-6.5));
    Text = 288.15 - 6.5*z;
    rho = Pext/(r*Text);

elseif z>=11 && z<20
    Pext = 22632.06*exp(-34.1632*(z-11) / 216.65);
    Text = 216.65;
    rho = Pext/(r*Text);

elseif z>=20 && z<32
    Pext = 5474.889*(216.65 / (216.65+(z-20)))^(34.1632);
    Text = 196.65 + z;
    rho = Pext/(r*Text);

elseif z>=32 && z<47
    Pext = 868.0187*(228.65 / (228.65 + 2.8*(z-32)))^(34.1632/2.8);
    Text = 139.05 + 2.8*z;
    rho = Pext/(r*Text);

elseif z>=47 && z<51
    Pext = 110.9063*exp(-34.1632*(z-47) / 270.65);
    Text = 270.65;
    rho = Pext/(r*Text);
    
elseif z>=51 && z<71
    Pext = 66.93887*(270.65 / (270.65-2.8*(z-51)))^(34.1632 / (-2.8));
    Text = 413.45 - 2.8*z;
    rho = Pext/(r*Text);
    
elseif z>=71 && z<84.852
    Pext = 3.956420*(214.65 / (214.65-2*(z-71)))^(34.1632 / (-2));
    Text = 356.65 - 2*z;
    rho = Pext/(r*Text);
    
elseif z>=84.852 && z<91
    Pext = exp( 0 * z^4 + 2.159582e-06 * z^3 - 4.836957e-04 * z^2 - 0.1425192 * z + 13.47530 );
    Text = 186.8673;
    rho = Pext/(r*Text);
    
elseif z>=91 && z<100
    Pext = exp( 0 * z^4 + 3.304895e-05 * z^3 - 0.009062730 * z^2 + 0.6516698 * z - 11.03037 );
    Text = 263.1905 - 76.3232 * sqrt(1 - ((z - 91) / (-19.9429))^2);
    rho = Pext/(r*Text);
    
elseif z>=100 && z<110
    Pext = exp( 0 * z^4 + 6.693926e-05	 * z^3 - 0.01945388 * z^2 + 1.719080 * z - 47.75030 );
    Text = 263.1905 - 76.3232 * sqrt(1 - ((z - 91) / (-19.9429))^2);
    rho = Pext/(r*Text);
    
elseif z>=110 && z<120
    Pext = exp( 0 * z^4 - 6.539316e-05	 * z^3 + 0.02485568 * z^2 - 3.223620 * z + 135.9355 );
    Text = 240 + 12 *(z - 110);
    rho = Pext/(r*Text);
    
elseif z>=120 && z<150
    Pext = exp(2.283506e-07 * z^4 - 1.343221e-04 * z^3 + 0.02999016 * z^2 - 3.055446 * z + 113.5764);
    Text = 1000 - 640 * exp(-0.01875*(z-120)*(6356.766+ 120)/(6356.766 +z));
    rho = Pext/(r*Text);
    
elseif z>=150 && z<200
    Pext = exp(1.209434e-08 * z^4 - 9.692458e-06 * z^3 + 0.003002041 * z^2 - 0.4523015 * z + 19.19151);
    Text = 1000 - 640 * exp(-0.01875*(z-120)*(6356.766+ 120)/(6356.766 +z));
    rho = Pext/(r*Text);
    
elseif z>=200 && z<300
    Pext = exp(8.113942e-10 * z^4 - 9.822568e-07 * z^3 + 4.687616e-4 * z^2 - 0.1231710 * z + 3.067409);
    Text = 1000 - 640 * exp(-0.01875*(z-120)*(6356.766+ 120)/(6356.766 +z));
    rho = Pext/(r*Text);
    
elseif z>=300 && z<=500
    Pext = exp(9.814674e-11 * z^4 - 1.654439e-7 * z^3 + 1.148115e-4 * z^2 - 0.05431334 * z - 2.011365);
    Text = 1000 - 640 * exp(-0.01875*(z-120)*(6356.766+ 120)/(6356.766 +z));
    rho = Pext/(r*Text);
    
else
    fprintf('Altitude out of range')
    disp(z);
end