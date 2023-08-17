function [A,B] = quad1D_linearization(x,u,params)
%QUAD2D_LINEARIZATION Summary of this function goes here
%   Detailed explanation goes here
gam = params.gam;
T = params.T;
tau = params.tau;

A = zeros(3);
B = zeros(3,1);

A(1,2) = 1;
A(2,2) = -gam;
A(2,3) = T*cos(x(3));
A(3,3) = -1/tau;

B(3) = 1/tau;
end

