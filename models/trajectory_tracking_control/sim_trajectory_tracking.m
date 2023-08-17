% Trajectory tracking with the 1D quadrotor. The frequency is adapted to
% the uncertainty of the model, which is learned from recorded
% trajectories.

clear
close all
rng('default')

%% Initialization
% Define quadrotor system
system = 'quad1D';
define_system

% Define simulation
dt_sim = 0.05;
Tsim = 6.3;
steps_per_sample = 10;
n_iter = 3;

% Define reference trajectory
tspan = 0:dt_sim:Tsim;
[x_ref,u_ref] = quad1D_reference_trajectory(tspan,params);
n_discretization_steps = 10;

% Training data
z_obs = [];
y_obs = [];

n_initial_samples = 20;
% Initial data from around the equilibrium
[z_obs,y_obs] = quad1D_samples(n_initial_samples,params);

% Initial data from around the reference trajectory
% n_initial_samples = n_discretization_steps;
% t_tmp = linspace(Tsim/n_initial_samples,Tsim,n_initial_samples);
% z_obs = zeros(n_initial_samples,n+m);
% y_obs = zeros(n_initial_samples,n);
% dev_from_ref_x = [0.5,0.5,0.2];     dev_from_ref_u = 0.2; 
% for i = 1:n_initial_samples
%     for j = 1:n
%         z_obs(i,j) = interp1(tspan,x_ref(j,:),t_tmp(i)) + rand*dev_from_ref_x(j);
%     end
%     for j = 1:m
%         z_obs(i,j+n) = interp1(tspan,u_ref(j,:),t_tmp(i)) + rand*dev_from_ref_u(j);
%     end
%     y_obs(i,:) = f(z_obs(i,1:n),z_obs(i,n+1:end),params);
% end

% Add noise
for i = 1:n
    y_obs(:,i) = y_obs(:,i) + sigma_n(i)*randn(n_initial_samples,1);
end

% Store trajectories
x_real = zeros(size(x_ref,1),size(x_ref,2),n_iter);
u_real = zeros(size(u_ref,1),size(u_ref,2),n_iter);

% Visualize initial training data
figure
scatter3(z_obs(:,1),z_obs(:,2),z_obs(:,3)); hold on
plot3(x_ref(1,:),x_ref(2,:),x_ref(3,:))
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');
figure
scatter3(z_obs(:,1),z_obs(:,3),z_obs(:,4)); hold on
plot3(x_ref(1,:),x_ref(3,:),u_ref(1,:))
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_3$','Interpreter','latex');
zlabel('$u$','Interpreter','latex');

Ts_vec = [0.25,0.5,1];
for iter = 1:n_iter
    %% Learning
    % Compute sampling time and updated controller based on available data
    % Fit GP models
    gpr = cell(n,1);
    for i = 1:n
        % Take subset of training samples (growing dataset)
        gpr{i} = fitrgp(z_obs,y_obs(:,i),'KernelFunction',...
            'ardsquaredexponential','Sigma',sigma_n(i),'ConstantSigma',true);
    end

    % Get maximum sampling time
    n_ctrl_updates = ceil(Tsim/Ts);
    K_all_current_iter = zeros(n_ctrl_updates,m,n);
    Ts_max_all_current_iter = zeros(1,n_discretization_steps);
    finished = 0;
    k = 1;
    while ~finished
        % Get (x,u) to linearize about
        t = k*(Tsim/n_discretization_steps);
        x_lin = zeros(1,n);
        u_lin = zeros(1,m);
        for i = 1:n
            x_lin(i) = interp1(tspan,x_ref(i,:),t);
        end
        for i = 1:m
            u_lin(i) = interp1(tspan,u_ref(i,:),t);
        end

        % Get linearized dynamics
        [A,B,Au,Bu] = get_linearized_dynamics(gpr,sigma_n,y_obs,x_lin,u_lin);
        [H,E,F] = transform_uncertainty(Au,Bu);       

        % Maximize sampling time
        eps_vec = logspace(-3,3,20);
        [Ts_max, K, tmp] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec);
        Ts_max_all_current_iter(k) = Ts_max;

        % If Ts_max == 0, increase dataset and start again
        if iter == 1 && Ts_max == 0
            fprintf("%i datapoints, infeasible at t = %.2f \n", size(z_obs,1), t);
            disp(Ts_max_all_current_iter);  
            k = 0;          % Start again

            % Get samples and add noise
            [z_obs_tmp,y_obs_tmp] = quad1D_samples(5,params);
            for i = 1:n
                y_obs_tmp(:,i) = y_obs_tmp(:,i) + sigma_n(i)*randn(5,1);
            end
            z_obs = [z_obs;z_obs_tmp];
            y_obs = [y_obs;y_obs_tmp];

            % Retrain with increased dataset
            gpr = cell(n,1);
            for i = 1:n
                % Take subset of training samples (growing dataset)
                gpr{i} = fitrgp(z_obs,y_obs(:,i),'KernelFunction',...
                    'ardsquaredexponential','Sigma',sigma_n(i),'ConstantSigma',true);
            end
        end

        % Check if finished
        if k == n_discretization_steps
            finished = 1;
            fprintf("%Finished with %i datapoints \n", size(z_obs,1));
        end
        k = k+1;
    end
    Ts = 0.5 * min(Ts_max_all_current_iter); 
    disp(Ts_max_all_current_iter);  
    fprintf("The sampling time in iteration %i is %.3f s\n",iter,Ts);
    Ts = Ts_vec(iter);

    %% Simulation  
    % Simulate system
    x_sim = []; u_sim = [];
    tspan_odeint = linspace(0,Ts*ceil(Tsim/Ts),(steps_per_sample-1)*ceil(Tsim/Ts)+1);
    x_k = x_ref(:,1)';
    y_obs_prev = y_obs;
    for k = 1:ceil(Tsim/Ts)
        disp(k)
        % Calculate (x,u) about which to linearize
        x_lin_k = zeros(1,n);
        u_lin_k = zeros(1,m);
        if k*Ts < Tsim      % Check if still in tspan
            for i = 1:n
                x_lin_k(i) = interp1(tspan,x_ref(i,:),k*Ts);
            end
            for i = 1:m
                u_lin_k(i) = interp1(tspan,u_ref(i,:),k*Ts);
            end
        else                % Last linearization point: end of trajectory
            x_lin_k = x_ref(:,end)';        
            u_lin_k = u_ref(:,end)';
        end

        % Get linearized dynamics
        [A,B,Au,Bu] = get_linearized_dynamics(gpr,sigma_n,y_obs_prev,x_lin_k,u_lin_k);
        % [A,B] = quad1D_linearization(x_lin,u_lin,params);
        % Au = zeros(size(A));    Bu = zeros(size(B));
        [H,E,F] = transform_uncertainty(Au,Bu);
        
        % Controller design       
        eps_vec = logspace(-2,2,10);
        [eta, K, eta_vec] = min_J_norm_bounded(A,B,H,E,F,Ts,J1,J2,eps_vec);

        % Calculate control input
        u = u_lin_k + K*(x_k - x_lin_k)';
    
        [t,x_new] = ode45(@(t,x) f(x,u,params), linspace(0,Ts,steps_per_sample), x_k);
        if k > 1
            x_new = x_new(2:end,:);
        end
    
        x_sim = [x_sim; x_new];
        u_sim = [u_sim; repmat(u',size(x_new,1),1)];
       
        % Store training sample
        z_tmp = [x_k,u];
        y_tmp = f(x_k,u,params)';
        
        % Add noise
        for i = 1:n
            y_tmp(i) = y_tmp(i) + sigma_n(i)*randn;
        end
        z_obs = [z_obs;z_tmp];
        y_obs = [y_obs;y_tmp];

        % Update x_k
        x_k = x_new(end,:);
    end

    % Store trajectories
    for i = 1:n
        x_real(i,:,iter) = interp1(tspan_odeint,x_sim(:,i),tspan)';
    end
    for i = 1:m
        u_real(i,:,iter) = interp1(tspan_odeint,u_sim(:,i),tspan);
    end
    
end

%% Plot
% Visualize final training data
figure
scatter3(z_obs(:,1),z_obs(:,2),z_obs(:,3)); hold on
plot3(x_ref(1,:),x_ref(2,:),x_ref(3,:))
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');
figure
scatter3(z_obs(:,1),z_obs(:,3),z_obs(:,4)); hold on
plot3(x_ref(1,:),x_ref(3,:),u_ref(1,:))
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_3$','Interpreter','latex');
zlabel('$u$','Interpreter','latex');

% Plot trajectories
figure
plot(x_ref(1,:),x_ref(3,:),'k--'); hold on
legend_entries = cell(1,n_iter+1);
legend_entries{1} = "Reference";
for i = 1:n_iter
    % Plot trajectory
    plot(x_real(1,:,i),x_real(3,:,i));
    legend_entries{i+1} = "Iteration " + num2str(i);
end
legend(legend_entries)
