% Simulate a nonlinear system which is stabilized arount a steady-state
% with a sampled-data controller. The simulation is repeated for different
% sampling times and initial conditions, and the amount of stable
% trajectories is counted.

clear
close all
rng('default')

%% Specify system model and simulation parameters
system = 'quad2D';
switch system
    case 'quad1D'
        n = 2; m = 1;
        f = quad1D_dynamics;
        x_s = [1,0,0];
        u_s = 0;        % Check again!
    case 'quad2D'
        n = 6; m = 2;
        f = @quad2D_dynamics;
        params = quad2D_params;
        x_s = [1,0,0,0,0,0];
        u_s = [0.1,0.1] * 9.81;
        K = [0.0022    0.0050   -0.0090   -0.0430   -0.0459   -0.0190;
            -0.0020   -0.0039   -0.0090   -0.0430    0.0442    0.0168];
        % [A,B] = quad2D_linearization(x_s,u_s,params);
        % K = -lqr(A,B,eye(6),eye(2));
    otherwise
        error('Wrong system type!')
end

% Specify sampling time and simulation time
Ts_vec = 0.05:0.05:0.5;
Tsim = 10;
steps_per_sample = 11;

% Specify how many initial conditions to try out
n_trials = 100;
x0_limits = [0.2,0,0.2,0,deg2rad(10),0];
stable_mat = ones(length(Ts_vec),n_trials);

%% Simulate system
for i = 1:length(Ts_vec)
    for j = 1:n_trials
        % Sample random initial state
        x_k = sample_x0(x_s,x0_limits);

        % Simulate system
        x_sim = []; u_sim = [];
        for k = 1:ceil(Tsim/Ts_vec(i))
            % Calculate input, which is kept constant
            u = K * (x_k - x_s)' + u_s';

            [t,x_new] = ode45(@(t,x) f(x,u,params), linspace(0,Ts_vec(i),steps_per_sample), x_k);
            if k > 1
                x_new = x_new(2:end,:);
            end
            x_sim = [x_sim; x_new];
            u_sim = [u_sim; repmat(u',size(x_new,1),1)];

            % Update x_k
            x_k = x_new(end,:);

            % Check if system has become unstable. If so, terminate simulation
            if norm([x_new(end,1), x_new(end,3)] - [x_s(1), x_s(3)]) > 1
                stable_mat(i,j) = 0;
                break
            end
        end
    end
end

%% Plot
% Stability as a function of sampling time
plot(Ts_vec,mean(stable_mat,2));

% Plot last state and input trajectories
if stable_mat(end,end)
    tspan = linspace(0,Ts_vec(end)*ceil(Tsim/Ts_vec(end)),(steps_per_sample-1)*ceil(Tsim/Ts_vec(end))+1);
    figure()
    plot(tspan,x_sim(:,1)); hold on
    plot(tspan,x_sim(:,3));
    plot(tspan,x_sim(:,5));
    legend('$x$','$z$','$\theta$','interpreter','latex')

    figure();
    plot(tspan,u_sim(:,1)); hold on
    plot(tspan,u_sim(:,2));
end

%% Helper functions
function x0_rand = sample_x0(x_s,x0_limits)
% Sample a random initial state. The i-th component is contained in 
% x0_rand(i) \in [x_s(i) - x0_limits(i), x_s(i) + x0_limits(i)]
    x0_rand = x_s;
    for i = 1:length(x_s)
        x0_rand(i) = x0_rand(i) + x0_limits(i) * (-1 + 2*rand);
    end
end
