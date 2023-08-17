function dx = quad2D_dynamics(x,u,params)
% 1D quadrotor model from Greeff and Schoellig, 2021.
% Define parameters
d = params.d;
m = params.m;
Iyy = params.Iyy;
g = 9.81;

if length(x) ~= 6
    error('State must have dimension of three!')
end

dx = zeros(6,1);
dx(1) = x(2);                               % x-velocity
dx(2) = -1/m * sin(x(5)*(u(1)+u(2)));       % x-acceleration
dx(3) = x(4);                               % z-velocity
dx(4) = 1/m * cos(x(5))*(u(1)+u(2)) - g;    % z-acceleration
dx(5) = x(6);                               % Angular velocity
dx(6) = 1/Iyy * d*(u(1)-u(2));              % Angular acceleration

end

