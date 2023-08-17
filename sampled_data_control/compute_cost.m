function J = compute_cost(x,u,x_s,u_s,J1,J2,Tsim)
% Get dimensions
[n_timesteps,n] = size(x);

% Calculate dt
dt = Tsim/n_timesteps;

J = 0;
for i = 1:n_timesteps
    J = J + ((x(i,:) - x_s) * J1 * (x(i,:) - x_s)' + (u(i,:) - u_s) * J2 * (u(i,:) - u_s)') * dt;
end

end

