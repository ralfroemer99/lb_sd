function [A,B] = quad2D_linearization(x,u,params)
%QUAD2D_LINEARIZATION Summary of this function goes here
%   Detailed explanation goes here
m = params.m;
d = params.d;
Iyy = params.Iyy;

A = zeros(6);
B = zeros(6,2);

A(1,2) = 1;
A(2,5) = -1/m * cos(x(5)) * (u(1)+u(2));
A(3,4) = 1;
A(4,5) = -1/m * sin(x(5)) * (u(1)+u(2));
A(5,6) = 1;

B(2,1) = -1/m * sin(x(5)); 
B(2,2) = B(2,1);
B(4,1) = 1/m * cos(x(5));
B(4,2) = B(4,1);
B(6,1) = d/Iyy;
B(6,2) = -d/Iyy;
end

