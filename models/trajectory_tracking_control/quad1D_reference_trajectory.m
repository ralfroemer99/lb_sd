function [x_ref, u_ref] = quad1D_reference_trajectory(tspan,p)
% Quad1D dynamics from Greeff and Schoellig, 2021
N = length(tspan);
x_ref = zeros(3,N);
u_ref = zeros(1,N);

% Define reference for position: p(t)=sin(t)
% k1 = 1;
% k2 = 1;
x1_ref = sin(tspan);            % x1
x2_ref = cos(tspan);            % x2
dx2_ref = -sin(tspan);          % dx2
ddx2_ref = -cos(tspan);         % ddx2
x3_ref = asin(1/p.T*(dx2_ref+p.gam*x2_ref));
dx3_ref = (1./sqrt(1-(1/p.T*(dx2_ref+p.gam*x2_ref)).^2)).*(1/p.T * (ddx2_ref + p.gam * dx2_ref));

x_ref(1,:) = x1_ref;
x_ref(2,:) = x2_ref;
x_ref(3,:) = asin(1/p.T*(dx2_ref+p.gam*x2_ref));

u_ref = x3_ref + p.tau*dx3_ref;

end

