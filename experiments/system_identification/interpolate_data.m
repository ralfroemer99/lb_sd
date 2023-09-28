function [data] = interpolate_data(data, dt)
    timestamps = data.iddata.timestamps;
    X = data.iddata.X;
    U = data.iddata.U;
    Xd = data.iddata.Xd;

    state_dim = size(X, 1);
    input_dim = size(U, 1);
    
    timestamps_interp = min(timestamps): dt: max(timestamps);
    X_interp = interp1(timestamps, X', timestamps_interp)';
    U_interp = interp1(timestamps, U', timestamps_interp)';
    Xd_interp = interp1(timestamps, Xd', timestamps_interp)';
    
    X_interp = reshape(X_interp, state_dim, []);
    U_interp = reshape(U_interp, input_dim, []);
    Xd_interp = reshape(Xd_interp, input_dim, []);
    
    % figure(1);
    % for i = 1:1:state_dim
    %     clf;
    %     hold on;
    %     plot(timestamps, states(i,:), '.', 'color', [0.5, 0.5, 0.5]);
    %     plot(timestamps_interp, states_interp(i,:), 'b-o');
    %     fprintf("Check inpterpolation. Press any key to continue.\n");
    %     legend('original data', 'interpolated data');
    %     xlabel('Time [s]');
    %     ylabel(sprintf('State %d', i));
    %     set(gca, 'fontsize', 12);
    %     waitforbuttonpress;
    % end
    % 
    % for i = 1:1:input_dim
    %     clf;
    %     hold on;
    %     plot(timestamps, inputs(i,:), '.', 'color', [0.5, 0.5, 0.5]);
    %     plot(timestamps_interp, inputs_interp(i,:), '-o', 'color', [0, 0, 1, 0.8]);
    %     fprintf("Check inpterpolation. Press any key to continue.\n");
    %     legend('original data', 'interpolated data');
    %     xlabel('Time [s]');
    %     ylabel(sprintf('Input %d', i));
    %     set(gca, 'fontsize', 12);
    %     waitforbuttonpress;
    % end
    
	data.iddata_interp.timestamps = timestamps_interp;
    data.iddata_interp.U = U_interp;
    data.iddata_interp.X = X_interp;
    data.iddata_interp.Xd = Xd_interp;
    data.iddata_interp.dt = dt;
    
    close all;
end