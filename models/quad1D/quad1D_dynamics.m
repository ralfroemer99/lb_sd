function dx = quad1D_dynamics(x,u,params)
% 1D quadrotor model from Greeff and Schoellig, 2021.
% Define parameters
tau = params.tau;
gam = params.gam;
T = params.T;

if length(x) ~= 3
    error('State must have dimension of three!')
end

dx = zeros(3,1);
dx(1) = x(2);                       % Linear velocity
dx(2) = T*sin(x(3)) - gam*x(2);     % Linear acceleration   
dx(3) = 1/tau * (u-x(3));           % Angular velocity
end

