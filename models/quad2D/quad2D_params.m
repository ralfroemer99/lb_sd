function p = quad2D_params()
%QUAD2D_PARAMS Summary of this function goes here
%   Detailed explanation goes here
p.d = 0.1;
p.m = 0.2;
p.Iyy = 1/12 * p.m * (2*p.d)^2;
end

