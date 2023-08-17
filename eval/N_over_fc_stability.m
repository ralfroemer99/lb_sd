% Learn the nonlinear dynamics from measured data and design a sampled-data
% controller for the linearized system.

% addpath '/home/ralf/mosek/10.0/toolbox/r2017aom'          % Ubuntu
% addpath 'C:\Program Files\Mosek\10.0\toolbox\r2017aom'    % Windows

close all
clear
rng('default')
warning('off','all')

save_path = "./data/N_over_fc_stability/sigman_0p1/fc_from_21_to_30";

f_vec = 21:1:30;
n_trials = 10;

%% Computation
% System model
system = 'quad2D';
define_system

N_vec = 100:25:1000; 

N_min_vec = zeros(length(f_vec),n_trials);
for k = 1:n_trials
    %% Learning
    fprintf("Seed %i \n", k)
    % Create noise-free samples
    switch system
        case 'quad1D'
            [z_obs,y_obs] = quad1D_samples(N_vec(end),params);
        case 'quad2D'
            [z_obs,y_obs] = quad2D_samples(N_vec(end),params);
    end

    % Add noise
    for i = 1:n
        y_obs(:,i) = y_obs(:,i) + sigma_n(i)*randn(N_vec(end),1);
    end
    
    % Train GP models with different amounts of data
    gpr = cell(1,length(N_vec));
    for i = 1:length(N_vec)
        gpr{i} = cell(1,length(n));
        for j = 1:n
            % Take subset of training samples (growing dataset)
            gpr{i}{j} = fitrgp(z_obs(1:N_vec(i),:),y_obs(1:N_vec(i),j),'KernelFunction',...
                'ardsquaredexponential','Sigma',sigma_n(j),'ConstantSigma',true);
        end
    end
    fprintf("Finished GP training \n")
    
    %% Optimization
    % Iterate over frequencies
    for i = 1:length(f_vec)
        % Iterate over dataset sizes
        for j = 1:length(N_vec)
            % Get linearized dynamics
            [A,B,Au,Bu] = get_linearized_dynamics(gpr{j},sigma_n,y_obs(1:N_vec(j),:),x_s,u_s);
            [H,E,F] = transform_uncertainty(Au,Bu);
            
            % Try to compute maximum sampling time
            eps_vec = logspace(-2,2,10);
            [Ts_max, K, Ts_vec] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec);

            % Check if sampling time is large enough
            if Ts_max >= 1/f_vec(i)                                 
                N_min_vec(i,k) = N_vec(j);
                fprintf("For f = %.2f, N_min = %i \n", f_vec(i), N_vec(j));
                break;   
            elseif j == length(N_vec)       % Reached end of dataset
                fprintf("For f = %.2f, N_min > %i \n", f_vec(i), N_vec(j));
            end
        end
    end
end

%% Plot
figure
N_min_mean = zeros(1,length(f_vec));
N_min_std = zeros(1,length(f_vec));
success_rate = zeros(1,length(f_vec));
for i = 1:length(f_vec)
    tmp = N_min_vec(i,:);
    N_min_mean(i) = mean(tmp(tmp ~= 0));
    N_min_std(i) = std(tmp(tmp ~= 0));
    success_rate(i) = length(tmp(tmp ~= 0)) / n_trials;
end
errorbar(f_vec,N_min_mean,N_min_std); hold on
figure
plot(f_vec, success_rate);

%% Save data
save(save_path + "f_vec",'f_vec');
save(save_path + "N_min_vec",'N_min_vec');
