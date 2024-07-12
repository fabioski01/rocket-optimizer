function Z = Fsyst2(y)
% Function that evaluates the two differential equations

Z = zeros(3,1);

Rt = 6371e3; % [m]
g = gravity(y(1));

Z(1) = y(2)*sin(y(3));

Z(2) = -g*sin(y(3));

Z(3) = -(g/y(2) - y(2)/(Rt+y(1)))*cos(y(3));

end