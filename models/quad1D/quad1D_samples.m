function [z_obs,y_obs] = quad1D_samples(N,params)
% 1D quadrotor model from Greeff and Schoellig, 2021.

% State limits
x_lim = zeros(3,2);
x_lim(1,:) = [-1.5,1.5];
x_lim(2,:) = [-1.5,1.5];
x_lim(3,:) = [-0.75,0.75];

% Input limits
u_lim = [-0.75,0.75];

% Create training inputs
x_obs = zeros(N,3);
u_obs = zeros(N,1);

x_obs(:,1) = x_lim(1,1) + (x_lim(1,2)-x_lim(1,1))*rand(N,1);
x_obs(:,2) = x_lim(2,1) + (x_lim(2,2)-x_lim(2,1))*rand(N,1);
x_obs(:,3) = x_lim(3,1) + (x_lim(3,2)-x_lim(3,1))*rand(N,1);
u_obs(:,1) = u_lim(1) + (u_lim(2)-u_lim(1))*rand(N,1);

z_obs = [x_obs,u_obs];

% Create training targets
y_obs = zeros(N,3);

for i = 1:N
    y_obs(i,:) = quad1D_dynamics(x_obs(i,:),u_obs(i),params);
end

end

