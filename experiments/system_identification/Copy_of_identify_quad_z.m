%%
clear
close all
clc
colors = matlab_plot_colors;
% rng('default')
rng(1,'twister')

%% Set Parameters
use_acc = true;
smooth_inputs = false;
control_freq = 30;    % Control frequency for interpolation
delta_t = 1 / control_freq;

% Specify sampled-data controller to compute
compute_Ts_max = false;
optimize_controller = false;
fc_des = 10; 
Q = eye(2);     R = 100;

N_how_many = 901;    % Set to [] to use all available data points
use_filtered_zdd_thrust = 1;

% Training data:
% file_ids = [3, 4, 5];   % 3
file_ids = [5];
file_path = "./experiments/data/training/track_sine_fc30trial";
% start_idx = [290 + 1, 290 + 1, 290 + 1];
% end_idx = [590 + 1, 1040 + 1, 1190 + 1]; 
start_idx = [290 + 1];
end_idx = [1190 + 1];  

% Trajectories:
% file_path = "./experiments/data/trajectories/setpoint_tracking_data_fc30";

X = [];
Y = [];

% selected_id = 2;
for k = 1:length(file_ids)
    %% Load data
    file = file_path + int2str(file_ids(k));
    % file = file_path;
    dronetrainingdata = load_drone_data(file, 1);

    % Extract relevant part of the trajectory
    dronetrainingdata = dronetrainingdata(start_idx(k):end_idx(k), :);
    
    %% Get data in the correct format
    time_id = 1;
    z_id = 2;
    zd_id = 3;
    zdd_id = 4;
    thrust_id = 5;
    thrust_delta_id = 6;
    
    u_id = thrust_id;
    
    % Set the input bounds
    u_min = -0.22;
    u_max = 0.25;
       
    iddata.timestamps = dronetrainingdata(:, time_id)';
    iddata.timestamps = iddata.timestamps - iddata.timestamps(1);
    x_ids = [z_id, zd_id];
    xd_ids = zdd_id;
    
    iddata.X = dronetrainingdata(:, x_ids)';
    iddata.U = dronetrainingdata(:, u_id)';
    iddata.Xd = dronetrainingdata(:, xd_ids)';
    
    iddata.U = invert_raw_to_delta(iddata.U);
    
    if u_id == thrust_id || u_id == thrust_delta_id
        iddata.U(iddata.U < u_min) = u_min;
        iddata.U(iddata.U > u_max) = u_max;
    end
    
    data.iddata = iddata;
    
    % interpolate data
    data = interpolate_data(data, delta_t);
    data = data.iddata_interp;
    data.dt = delta_t;

    % Determine z-acc by numerical differentiation
    zdd = diff(data.X(2, :) / delta_t);
    zdd = [zdd, zdd(end)];

    % Mean filter acceleration and thrust
    window_size = 5;
    zdd_filtered = medfilt1(zdd,window_size);
    zdd_filtered = [zdd_filtered(2), zdd_filtered(2:end-1), zdd_filtered(end-1)];
    data.U_filtered = medfilt1(data.U,window_size);
    data.U_filtered = [data.U_filtered(2), data.U_filtered(2:end-1), data.U_filtered(end-1)];

    % iddata.u = smoothdata(dronetrainingdata(:, input_id)');
    
    %% Plot with integration check
    integrate_vel = cumsum(data.X(2, :) * delta_t) + data.X(1, 1);
    integrate_acc = cumsum(data.Xd * delta_t) + data.X(2, 1);
    diff_pos = diff(data.X(1, :)) / delta_t;
    
    % Position
    subplot(4, 1, 1)
    hold on
    plot(data.timestamps, data.X(1, :), "Color", colors(1,:), "DisplayName", "Meas. $z$")
    plot(data.timestamps(1:end), integrate_vel, "Color", colors(2,:), "DisplayName", "Meas. $\dot{z}$ integrated")
    legend('Interpreter','latex')
    xlabel("$t$ in s", 'Interpreter','latex')
    ylabel("$z$ in m", 'Interpreter','latex')
    title("Position")
    
    % Velocity
    subplot(4, 1, 2)
    hold on
    plot(data.timestamps, data.X(2, :), "Color", colors(1,:), "DisplayName", "Meas. $\dot{z}$")
    plot(data.timestamps(1:end), integrate_acc, "Color", colors(2,:), "DisplayName", "Meas. $\ddot{z}$ integrated")
    plot(data.timestamps(2:end), diff_pos, "Color", colors(3,:), "DisplayName", "Meas. $z$ differentiated")
    legend('Interpreter','latex')
    xlabel("$t$ in s", 'Interpreter','latex')
    ylabel("$\dot{z}$ in $\frac{\mathrm{m}}{\mathrm{s}}$", 'Interpreter','latex')
    title("Velocity")
    
    % Measured Acceleration vs. Differentiated Position
    subplot(4, 1, 3)
    hold on
    % plot(data.timestamps, data.Xd, "Color", colors(1,:), "DisplayName", "Meas. $\ddot{z}$")
    plot(data.timestamps, zdd, "Color", colors(2,:), "DisplayName", "Meas. $\dot{z}$ diff.")
    plot(data.timestamps, zdd_filtered, "Color", colors(3,:), "DisplayName", "Meas. $\dot{z}$ diff. filtered")
    legend('Interpreter','latex')
    xlabel("$t$ in s", 'Interpreter','latex')
    ylabel("$\ddot{z}$ in $\frac{\mathrm{m}}{\mathrm{s}^2}$", 'Interpreter','latex')
    title("Acceleration")
    hold off

    % Thrust
    subplot(4, 1, 4)
    hold on
    plot(data.timestamps, data.U, "Color", colors(1,:), "DisplayName", "Meas. $T(t)$")
    plot(data.timestamps, data.U_filtered, "Color", colors(2,:), "DisplayName", "Meas. $T(t)$ filtered")
    legend('Interpreter','latex')
    xlabel("$t$ in s", 'Interpreter','latex')
    ylabel("$T(t)$ in $\mathrm{N}$", 'Interpreter','latex')
    title("Thrust")
    hold off
    
    if use_filtered_zdd_thrust   
        Y = [Y; 
             zdd_filtered'];
        X = [X;
             data.X(1, :)', data.X(2, :)', data.U_filtered'];
    else
        Y = [Y; 
             zdd'];
        X = [X;
             data.X(1, :)', data.X(2, :)', data.U'];
    end
end

%% Identify model iteratively using BLR
N = size(X,1);
% sigma_n = [1,1,0.1];
% Calculate meaningful noise variance: From std(Y' - data.inputs / 0.033)
sigma_n = std(Y' - data.U / 0.033)^2;
% sigma_n = 1;
mu_theta_vec = [0; 0; 1/0.033];
Sigma_theta_cell = cell(1,N+1);
Sigma_theta_cell{1}= eye(3);
Sigma_theta_vec = zeros(3,N+1);
Sigma_theta_vec(:,1) = [1; 1; 100];

reshuffled_indices = randperm(N,N);
for i = 1:N
    % Extract a random datapoint
    j = reshuffled_indices(i);
    % Update mean vector and covariance matrix
    Sigma_theta_cell{i+1} = (Sigma_theta_cell{i}^(-1) + diag(sigma_n.^(-2)) * X(j,:)'*X(j,:))^(-1);
    mu_theta_vec(:,i+1) = Sigma_theta_cell{i+1}*(Sigma_theta_cell{i}\mu_theta_vec(:,i) + diag(sigma_n.^(-2)) * X(j,:)' * Y(j));
    Sigma_theta_vec(:,i+1) = diag(Sigma_theta_cell{i+1});
end

% Get parameter estimate
C = [0 1 0; 0 0 0];
C(2,:) = mu_theta_vec(:,N_how_many+1);
A = C(:,1:2)       
B = C(:,3)

% Convert paramter distribution to confidence intervals and reparemeterize
C_hat = zeros(2,3);
C_hat(2,:) = sqrt(chi2inv(0.95,3)) * sqrt(diag(Sigma_theta_cell{N_how_many+1}))';
A_hat = C_hat(:,1:2)   
B_hat = C_hat(:,3)
[H,E,F] = transform_uncertainty(A_hat,B_hat);

% Compute maximum sampling time
if compute_Ts_max
    eps_vec = logspace(-3,3,20);
    [Ts_max, K, ~] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec)
    if Ts_max ~= 0
        fc_min = 1/Ts_max
    end
end

% Optimize controller
if optimize_controller
    eps_vec = logspace(-3,3,20);
    [~, K_opt, eta_vec] = min_J_norm_bounded(A,B,H,E,F,1/fc_des,Q,R,eps_vec);
    K_opt
end

% Plot parameter mean
figure
for i = 1:3
    subplot(1,3,i)
    tmp = mu_theta_vec(i,:); plot(0:N,tmp(:)); hold on
    legend("$\mu_" + i + "$",'interpreter','latex')
end

% Plot parameter covariances
figure
for i = 1:3
    subplot(1,3,i)
    plot(0:N,Sigma_theta_vec(i,:)); hold on
    legend("$\Sigma_{" + i + i + "}$",'interpreter','latex')
end

% Plot parameter mean with \pm 2 std.
figure
for i = 1:3
    subplot(3,1,i)
    tmp_mean = mu_theta_vec(i,:);
    tmp_std = sqrt(Sigma_theta_vec(i,:));
    [tmp1, tmp2] = shaded_plot_mean_std(0:N, tmp_mean, tmp_std);
    fill(tmp1, tmp2, colors(i,:),'EdgeColor',colors(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.2,'HandleVisibility','off'); hold on
    plot(0:N,tmp_mean, 'Color', colors(i,:), 'LineWidth', 1);
end

% Fit integrator model manually
params = X \ Y;
A = [0.0, 1; params(1), params(2)];
B = [0.0; params(3)];
