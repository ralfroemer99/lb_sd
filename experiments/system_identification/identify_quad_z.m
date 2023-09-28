%%
clear
close all
clc

%% Set Parameters
use_acc = true;
smooth_inputs = false;
control_freq = 30;    % Control frequency for interpolation
delta_t = 1 / control_freq;

compute_Ts_max = true;

file_ids = 4;   % 3
% file_path = "./experiments/data/trajectories/setpoint_tracking_data_fc30";
file_path = "./experiments/data/training/track_sine_fc30trial4";
start_idx = 290 + 1;    % 290 + 1
end_idx = 890 + 1;     % 590 + 1

X = [];
Y = [];

% selected_id = 2;
for selected_id = 1:length(file_ids)
    %% Load data
    % file = file_path + int2str(file_ids(selected_id));
    file = file_path;
    dronetrainingdata = load_drone_data(file, 1);

    % Extract relevant part of the trajectory
    dronetrainingdata = dronetrainingdata(start_idx:end_idx, :);
    
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
      
    % iddata.u = smoothdata(dronetrainingdata(:, input_id)');
    
    %% Plot with integration check
    integrate_vel = cumsum(data.X(2, :) * delta_t) + data.X(1, 1);
    integrate_acc = cumsum(data.Xd * delta_t) + data.X(2, 1);
    diff_pos = diff(data.X(1, :)) / delta_t;
    diff_vel = diff(data.X(2, :)) / delta_t;
    
    % Position
    subplot(3, 1, 1)
    hold on
    plot(data.timestamps, data.X(1, :), "Color", "r", "DisplayName", "Meas. $z$")
    plot(data.timestamps(1:end), integrate_vel, "Color", "b", "DisplayName", "Meas. $\dot{z}$ integrated")
    legend('Interpreter','latex')
    xlabel("$t$ in s", 'Interpreter','latex')
    ylabel("$z$ in m", 'Interpreter','latex')
    title("Position")
    
    % Velocity
    subplot(3, 1, 2)
    hold on
    plot(data.timestamps, data.X(2, :), "Color", "r", "DisplayName", "Meas. $\dot{z}$")
    plot(data.timestamps(1:end), integrate_acc, "Color", "b", "DisplayName", "Meas. $\ddot{z}$ integrated")
    plot(data.timestamps(2:end), diff_pos, "Color", "g", "DisplayName", "Meas. $z$ differentiated")
    legend('Interpreter','latex')
    xlabel("$t$ in s", 'Interpreter','latex')
    ylabel("$\dot{z}$ in $\frac{\mathrm{m}}{\mathrm{s}}$", 'Interpreter','latex')
    title("Velocity")
    
    % Measured Acceleration vs. Differentiated Position
    subplot(3, 1, 3)
    hold on
    plot(data.timestamps, data.Xd, "Color", "r", "DisplayName", "Meas. $\ddot{z}$")
    plot(data.timestamps(2:end), diff_vel, "Color", "g", "DisplayName", "Meas. $\dot{z}$ differentiated")
    legend('Interpreter','latex')
    xlabel("$t$ in s", 'Interpreter','latex')
    ylabel("$\ddot{z}$ in $\frac{\mathrm{m}}{\mathrm{s}^2}$", 'Interpreter','latex')
    title("Acceleration")
    hold off
    
    Y = [Y; 
         zdd'];
    X = [X;
         data.X(1, :)', data.X(2, :)', data.U'];
end

%% Identify model iteratively using BLR
N = size(X,1);
% sigma_n = [1,1,0.1];
% Calculate meaningful noise variance: From std(Y' - data.inputs / 0.033)
sigma_n = std(Y' - data.U / 0.033)^2;
mu_theta_vec = [0; 0; 1/0.033];
Sigma_theta_cell = cell(1,N+1);
Sigma_theta_cell{1}= eye(3);
Sigma_theta_vec = zeros(3,N+1);
Sigma_theta_vec(:,1) = [1; 1; 10];
for i = 1:N
    % Update mean vector and covariance matrix
    Sigma_theta_cell{i+1} = (Sigma_theta_cell{i}^(-1) + diag(sigma_n.^(-2)) * X(i,:)'*X(i,:))^(-1);
    mu_theta_vec(:,i+1) = Sigma_theta_cell{i+1}*(Sigma_theta_cell{i}\mu_theta_vec(:,i) + diag(sigma_n.^(-2)) * X(i,:)' * Y(i));
    Sigma_theta_vec(:,i+1) = diag(Sigma_theta_cell{i+1});
end

% Get parameter estimate
C = [0 1 0; 0 0 0];
C(2,:) = mu_theta_vec(:,end);
A = C(:,1:2)       
B = C(:,3)

% Convert paramter distribution to confidence intervals and reparemeterize
C_hat = zeros(2,3);
C_hat(2,:) = sqrt(chi2inv(0.95,3)) * sqrt(diag(Sigma_theta_cell{end}))';
A_hat = C_hat(:,1:2)   
B_hat = C_hat(:,3)
[H,E,F] = transform_uncertainty(A_hat,B_hat);

% Compute maximum sampling time
if compute_Ts_max
    eps_vec = logspace(-3,3,20);
    [Ts_max, K, Ts_vec] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec);
    Ts_max
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

% Fit integrator model manually
params = X \ Y;
A = [0.0, 1; params(1), params(2)];
B = [0.0; params(3)];

