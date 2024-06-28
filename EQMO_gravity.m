function g_h = gravity(h)
g0 = 9.81; % [m/s^2]
Rt = 6371e3; % [m]
g_h = g0/(1+h/Rt)^2;
end