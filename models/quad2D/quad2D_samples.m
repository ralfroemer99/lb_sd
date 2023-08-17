function [z_obs,y_obs] = quad2D_samples(N,params)
% 1D quadrotor model from Greeff and Schoellig, 2021.

% State limits
x_lim = zeros(6,2);
x_lim(1,:) = [0,2];         % x
x_lim(2,:) = [-5,5];        % dx
x_lim(3,:) = [0,2];         % z
x_lim(4,:) = [-5,5];        % dz
x_lim(5,:) = [-pi/2,pi/2];  % theta
x_lim(6,:) = [-5,5];        % dtheta

% Input limits
u_lim = zeros(2,2);
u_lim(1,:) = [0,2];
u_lim(2,:) = [0,2];

% Create training inputs
x_obs = zeros(N,6);
for i = 1:6
    x_obs(:,i) = x_lim(i,1) + (x_lim(i,2)-x_lim(i,1))*rand(N,1);
end

u_obs = zeros(N,2);
u_obs(:,1) = u_lim(1,1) + (u_lim(1,2)-u_lim(1,1))*rand(N,1);
u_obs(:,2) = u_obs(:,1) + 0.1*randn(N,1);       % Only small thrust differences

z_obs = [x_obs,u_obs];

% Create training targets
y_obs = zeros(N,6);

for i = 1:N
    y_obs(i,:) = quad2D_dynamics(x_obs(i,:),u_obs(i,:),params);
end

end

