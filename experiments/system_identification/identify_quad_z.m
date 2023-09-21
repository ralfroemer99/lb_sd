%%
clear
close all
clc

%% Set Parameters
catch_data = true;
use_acc = false;
smooth_inputs = false;
index = 1;              % Setting the index for shifting the timestamp data
control_freq = 30.0;    % Control frequency for interpolation
delta_t = 1 / control_freq;
delay = 0;
discrete = false;

file_path = "./experiments/data/from_lukas/trajectories/drone_training_data_catch_";
file_ids = 1:4;

X = [];
Y = [];

% selected_id = 2;
for selected_id = 1:length(file_ids)

    %% Load data
    file = file_path + int2str(file_ids(selected_id));
    dronetrainingdata = load_drone_data(file, catch_data);
    
    %% Get data in the correct format
    time_id = 1;
    z_id = 2;
    zd_id = 3;
    thrust_id = 4;
    thrust_delta_id = 5;
    z_ref_id = 6;
    if ~catch_data
        za_id = 7;
    else
        za_id = 6;
    end
    
    input_id = thrust_delta_id;
    input_name = "";
    switch input_id
        case thrust_id
            input_name = "thrust";
        case thrust_delta_id
            input_name = "delta\_thrust";
        case z_ref_id
            input_name = "z\ref";
        otherwise
            disp("Invalid input_id provided!")
            exit
    end
    
    % Set the input bounds
    if input_id == thrust_delta_id
        u_min = - 0.22;
        u_max = 0.25;
    elseif input_id == thrust_id
        u_min = 2e4;
        u_max = 2^16 - 1;
    end
    
    % incorrect data recorded. Fixing this here by recalculating delta thrust
    % from thrust
    if catch_data  && input_id == thrust_delta_id
        input_id = thrust_id;
    end
    
    iddata.timestamps = dronetrainingdata(:, time_id)';
    state_ids = [z_id, zd_id];
    if use_acc
        state_ids = [state_ids, za_id];
    end
    
    iddata.states = dronetrainingdata(:, state_ids)';
    iddata.inputs = dronetrainingdata(:, input_id)';
    
    if catch_data
        iddata.inputs = invert_raw_to_delta(iddata.inputs);
    end
    
    if input_id == thrust_id || input_id == thrust_delta_id
        iddata.inputs(iddata.inputs < u_min) = u_min;
        iddata.inputs(iddata.inputs > u_max) = u_max;
    end
    
    data.iddata = iddata;
    
    % interpolate data
    data = interpolate_data(data, delta_t);
    data = data.iddata_interp;
    data.dt = delta_t;
    
    % offset data by certain number of timesteps
    data.timestamps = data.timestamps(index:end);
    data.states = data.states(:, index:end);
    data.inputs = data.inputs(index:end);
    
    % iddata.inputs = smoothdata(dronetrainingdata(:, input_id)');
    
    %% Integation check
    if use_acc
        integrate_acc = cumsum(data.states(3, :) * delta_t) + data.states(2, 1);
    end
    integrate_vel = cumsum(data.states(2, :) * delta_t) + data.states(1, 1);
    diff_pos = diff(data.states(1, :)) / delta_t;
    diff_vel = diff(data.states(2, :)) / delta_t;
    
    subplot(2, 2, 1)
    hold on
    plot(data.timestamps, data.states(1, :), "Color", "r", "DisplayName", "z")
    plot(data.timestamps(1:end), integrate_vel, "Color", "b", "DisplayName", "z\_vel\_int")
    legend
    xlabel("t [s]")
    ylabel("z [m]")
    title("Integrate Velocity")
    hold off
    
    subplot(2, 2, 2)
    hold on
    plot(data.timestamps, data.states(2, :), "Color", "r", "DisplayName", "z\_dot")
    if use_acc
        plot(data.timestamps(1:end), integrate_acc, "Color", "b", "DisplayName", "z\_acc\_int")
    end
    legend
    xlabel("t [s]")
    ylabel("z\_vel [m/s]")
    title("Integrate Acc.")
    hold off
    
    subplot(2, 2, 3)
    hold on
    plot(data.timestamps(2:end), diff_vel, "Color", "b", "DisplayName", "z\_vel\_diff")
    if use_acc
        plot(data.timestamps, data.states(3, :), "Color", "r", "DisplayName", "z\_acc")
    end
    legend
    xlabel("t [s]")
    ylabel("z\_acc [m/s^2]")
    title("Differentiate Velocity")
    hold off
    
    subplot(2, 2, 4)
    hold on
    plot(data.timestamps(2:end), diff_pos, "Color", "b", "DisplayName", "z\_diff")
    plot(data.timestamps, data.states(2, :), "Color", "r", "DisplayName", "z\_vel")
    legend
    xlabel("t [s]")
    ylabel("z\_vel [m/s]")
    title("Differentiate Position")
    hold off
    
    %% Fit LTI system
    % idmodel = identify_ss_model(data, delay, discrete, use_acc)
    % idmodel_vel = identify_vel_dyn(data, delay, discrete)

    % determine z-acc
    acc = diff(data.states(2, :) / delta_t);
    acc = [acc, acc(end)];
    
    Y = [Y; 
         acc'];
    X = [X;
         data.states(1, :)', data.states(2, :)', data.inputs'];

end

%% Identify model iteratively using BLR
N = size(X,1);
% mu_theta = zeros(3,1);
% Sigma_theta = ones(3);
sigma_n = [1,1,0.1];
mu_theta_vec = zeros(3,N+1);
Sigma_theta_vec = cell(1,N);
Sigma_theta_vec{1}= eye(3);
for i = 1:N
    % Update mean vector and covariance matrix
    Sigma_theta_vec{i+1} = (Sigma_theta_vec{i}^(-1) + diag(sigma_n.^(-2)) * X(i,:)'*X(i,:))^(-1);
    mu_theta_vec(:,i+1) = Sigma_theta_vec{i+1}*(Sigma_theta_vec{i}\mu_theta_vec(:,i) + diag(sigma_n.^(-2)) * X(i,:)' * Y(i));
end

% Get paramter estimate
C = [0 1 0; 0 0 0];
C(2,:) = mu_theta_vec(:,end);
A = C(:,1:2);       B = C(:,3);

% Convert paramter distribution to confidence intervals and reparemeterize
C_hat = zeros(2,3);
C_hat(2,:) = sqrt(chi2inv(0.95,3)) * sqrt(diag(Sigma_theta_vec{end}))';
A_hat = C_hat(:,1:2);   B_hat = C_hat(:,3);
[H,E,F] = transform_uncertainty(A_hat,B_hat);

% Compute maximum sampling time
% eps_vec = [logspace(-3,3,10); logspace(-3,3,10)];
eps_vec = logspace(-3,3,10);
[Ts_max, K, Ts_vec] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec);

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
    tmp = Sigma_theta_vec(i,i,:); plot(0:N,tmp(:)); hold on
    legend("$\Sigma_{" + i + i + "}$",'interpreter','latex')
end

% Fit integrator model manually
params = X \ Y;
A = [0.0, 1; params(1), params(2)];
B = [0.0; params(3)];

