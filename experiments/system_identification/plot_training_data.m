%%
clear
close all
clc
colors = matlab_plot_colors;
rng('default')
% rng(1,'twister')

%% Set Parameters
control_freq = 30;    % Control frequency for interpolation
delta_t = 1 / control_freq;

use_filtered_zdd_thrust = 1;

% Training data:
% file_ids = [3, 4, 5];   % 3
file_ids = [3];
file_path = "./experiments/data/training/track_sine_fc30_different_amp_trial";
% start_idx = [290 + 1, 290 + 1, 290 + 1];
% end_idx = [590 + 1, 1040 + 1, 1190 + 1]; 
start_idx = [290 + 1];
end_idx = [1190 + 1];  

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
    z_des = 1 + linspace(0.1,0.5,length(data.timestamps)).*sin(linspace(0,12*pi,length(data.timestamps)));
    plot(data.timestamps, z_des, "Color", 'black', "DisplayName", "$z_\mathrm{des}$")
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
