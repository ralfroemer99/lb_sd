%%
clear
close all
clc
colors = matlab_plot_colors;
rng('default')
% rng(1,'twister')

control_freq = 30;    % Control frequency for interpolation
delta_t = 1 / control_freq;
m = 0.035;

% Specify sampled-data controller to compute
N_vec = [250, 500, 750];
% n_trials = 10;
% N_vec = 900;
n_trials = 1;

use_filtered = 1;

optimize_controller = 1;    
% Q = eye(2);     R = 1;
Q = diag([100, 1]);     R = 1;
fc_des = [60, 30, 15, 10];
% fc_des = [10];
K_opt_vec = zeros(2,length(N_vec),length(fc_des),n_trials);

% Training data:
file_ids = [3];
file_path = "./experiments/data/training/track_sine_fc30_different_amp_trial";
start_idx = [290 + 1];
end_idx = [1190 + 1];

[X,Y,zdd] = identify_quad_z(file_path,file_ids,start_idx,end_idx,delta_t,use_filtered);

%% Identify model iteratively using BLR
N = size(X,1);
% Calculate meaningful noise variance: From std(Y' - data.inputs / m)
sigma_n = 2 * std(Y - X(:,3) / m);

Ts_max_vec = zeros(n_trials,length(N_vec));

for j = 1:n_trials
    mu_theta_vec = [0; 0; 1/m];
    Sigma_theta_cell = cell(1,N+1);
    Sigma_theta_cell{1}= eye(3);
    Sigma_theta_vec = zeros(3,N+1);
    Sigma_theta_vec(:,1) = [1; 1; 100];

    reshuffled_indices = randperm(N,N);
    for i = 1:N
        % Extract a random datapoint
        i_rs = reshuffled_indices(i);
        % Update mean vector and covariance matrix
        Sigma_theta_cell{i+1} = (Sigma_theta_cell{i}^(-1) + diag(sigma_n.^(-2)) * X(i_rs,:)'*X(i_rs,:))^(-1);
        mu_theta_vec(:,i+1) = Sigma_theta_cell{i+1}*(Sigma_theta_cell{i}\mu_theta_vec(:,i) + diag(sigma_n.^(-2)) * X(i_rs,:)' * Y(i_rs));
        Sigma_theta_vec(:,i+1) = diag(Sigma_theta_cell{i+1});
    end

    for k = 1:length(N_vec)
        % Get parameter estimate
        C = [0 1 0; 0 0 0];
        C(2,:) = mu_theta_vec(:,N_vec(k)+1);
        A = C(:,1:2);
        B = C(:,3);

        % Convert paramter distribution to confidence intervals and reparemeterize
        C_hat = zeros(2,3);
        C_hat(2,:) = sqrt(chi2inv(0.99,3)) * sqrt(diag(Sigma_theta_cell{N_vec(k)+1}))';
        A_hat = C_hat(:,1:2);
        B_hat = C_hat(:,3);
        [H,E,F] = transform_uncertainty(A_hat,B_hat);

        % Compute maximum sampling time
        eps_vec = logspace(-3,3,20);
        [Ts_max, K, ~] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec);
        Ts_max_vec(j,k) = Ts_max;
        fprintf("Trial %0.f, N = %.0f, Ts_max = %.4f \n", j, N_vec(k), Ts_max)

        % Optimize controller
        if optimize_controller
            eps_vec = logspace(-3,3,20);
            for l = 1:length(fc_des)
                disp(l)
                [~, K_opt, eta_vec] = min_J_norm_bounded(A,B,H,E,F,1/fc_des(l),Q,R,eps_vec);
                K_opt_vec(:,k,l,j) = K_opt;
            end
        end
    end

    % Plot parameter mean
    figure(1)
    for i = 1:3
        subplot(3,1,i)
        tmp = mu_theta_vec(i,:); 
        plot(0:N,tmp(:), "DisplayName", "Trial = " + j); hold on
        legend()
    end

    % Plot parameter covariances
    figure(2)
    for i = 1:3
        subplot(3,1,i)
        plot(0:N,Sigma_theta_vec(i,:), "DisplayName", "Trial = " + j); hold on
        legend()
    end

    % Plot with error bars
    figure(3)
    for i = 1:3
        subplot(3,1,i)
        tmp_mean = mu_theta_vec(i,:);
        tmp_std = sqrt(Sigma_theta_vec(i,:));
        [tmp1, tmp2] = shaded_plot_mean_std(0:N, tmp_mean, tmp_std);
        fill(tmp1, tmp2, colors(i,:),'EdgeColor',colors(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.2,'HandleVisibility','off'); hold on
        plot(0:N,tmp_mean, 'Color', colors(i,:), 'LineWidth', 1);
    end
end

% Plot fc_min
if length(N_vec) > 1
    success_rate = zeros(1,length(N_vec));
    % Remove results where Ts_max = 0
    fc_min_mean = zeros(1,length(N_vec));
    fc_min_std = zeros(1,length(N_vec));
    for k = 1:length(N_vec)
        tmp1 = Ts_max_vec(:,k);
        tmp2 = tmp1(tmp1 >= 0.001);
        fc_min_mean(k) = mean(1./tmp2);
        fc_min_std(k) = std(1./tmp2);
        success_rate(k) = length(tmp2) / length(tmp1);
    end    
  
    % Plot minimum control frequency
    figure()
    errorbar(N_vec,fc_min_mean,fc_min_std);
    % plot(N_vec, fc_min_mean);
    xlabel('Number of data points $N$','Interpreter','latex');
    ylabel('Minimum control frequency $f_{c,min}$','Interpreter','latex');
    title('Minimum control frequency')
    xlim([N_vec(1) - 25, N_vec(end) + 25])

    % Plot success rate, i.e., how many times a stabilizing controller 
    % could be found
    figure()
    scatter(N_vec, success_rate,'r','Marker','o');
    xlabel('Number of data points $N$','Interpreter','latex');
    ylabel('Success rate','Interpreter','latex');
    title('Success rate in computing a robustly stabilizing controller')
    xlim([N_vec(1) - 25, N_vec(end) + 25])
    ylim([-0.05 1.05])
end
