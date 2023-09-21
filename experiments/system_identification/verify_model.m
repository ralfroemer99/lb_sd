function [] = verify_model(model, data, discrete, x_desired, K)
    % Get system matrices
    A = model.A;
    B = model.B;

    if nargin < 4
        follow_inputs = true;
    else
        follow_inputs = false;
        A_k = (A + B * K);
    end
    
    % Get data
    timestamps = data.timestamps;
    states = data.states;
    inputs = data.inputs;
    state_dim = size(data.states, 1)
    % inputs = 0.0 * inputs
    
    %% Simulate model using control input data
    X_sim = zeros(size(states));
    X_sim(:, 1) = states(:, 1);

    mean(diff(timestamps))

    if discrete    
        for i = 2 : length(timestamps)
            if follow_inputs
                X_sim(:, i) = A * X_sim(:, i - 1) + B * inputs(:, i - 1);
            else
                X_sim(:, i) = A_k * (x_desired - X_sim(:, i - 1));
                inputs(:, i) = K * (x_desired - X_sim(:, i - 1));
            end
        end
    else
        X_sim = lsim(model, inputs, timestamps, X_sim(:, 1))';
    end

    %% Compare simulator system and identified system
    figure;
    clf;
    if state_dim > 1
        subplot(state_dim,1,1);
        hold on;
        box on;
        plot(timestamps, states(1, :), 'r*');
        plot(timestamps, X_sim(1, :), 'b*-');
        % plot(timestamps, inputs(1, :), '--x', 'color', [1, 1, 1] * 0.5);
        ylabel('Position [m]');
	    set(gca, 'fontsize', 12);
        xticks({});
    end
    
    if state_dim > 1
        subplot(state_dim,1,2);
        hold on;
        box on;
        plot(timestamps, states(2, :), 'r*');
        plot(timestamps, X_sim(2, :), 'b*-');
        if state_dim < 3
            plot(timestamps, inputs(1, :), '--x', 'color', [1, 1, 1] * 0.5);
        end
        ylabel('Velocity [m/s]');
    else
        subplot(state_dim,1,1);
        hold on;
        box on;
        plot(timestamps, states(1, :), 'r*');
        plot(timestamps, X_sim(1, :), 'b*-');
        plot(timestamps, inputs(1, :), '--x', 'color', [1, 1, 1] * 0.5);
        ylabel('Velocity [m/s]');
    end
    
    if state_dim > 2
        subplot(state_dim,1,3);
        hold on;
        box on;
        plot(timestamps, states(3, :), 'r*');
        plot(timestamps, X_sim(3, :), 'b*-');
        plot(timestamps, inputs(1, :), '--x', 'color', [1, 1, 1] * 0.5);
        ylabel('Acc [m/s2]');
    end

    xlabel('Time [s]')
    legend('Actual Output', 'Model Output', 'Input', 'location', 'sw');
	set(gca, 'fontsize', 12);
    hold off
end