function [data] = interpolate_data(data, dt)
    timestamps = data.iddata.timestamps;
    states = data.iddata.states;
    inputs = data.iddata.inputs;

    state_dim = size(states, 1);
    input_dim = size(inputs, 1);
    
    timestamps_interp = min(timestamps): dt: max(timestamps);
    states_interp = interp1(timestamps, states', timestamps_interp)';
    inputs_interp = interp1(timestamps, inputs', timestamps_interp)';
    
    states_interp = reshape(states_interp, state_dim, []);
    inputs_interp = reshape(inputs_interp, input_dim, []);
    
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
    data.iddata_interp.inputs = inputs_interp;
    data.iddata_interp.states = states_interp;
    data.iddata_interp.dt = dt;
    
    close all;
end