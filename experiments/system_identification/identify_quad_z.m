function [X,Y,zdd] = identify_quad_z(file_path,file_ids,start_idx,end_idx,delta_t,use_filtered)

window_size = 5;

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
    zdd_filtered = medfilt1(zdd,window_size);
    zdd_filtered = [zdd_filtered(2), zdd_filtered(2:end-1), zdd_filtered(end-1)];
    data.U_filtered = medfilt1(data.U,window_size);
    data.U_filtered = [data.U_filtered(2), data.U_filtered(2:end-1), data.U_filtered(end-1)];

    if use_filtered   
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
