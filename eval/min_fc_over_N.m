% Learn the nonlinear dynamics from measured data and compute the minimum
% control frequency required to stabilize the uncertain linearized system.

% addpath '/home/ralf/mosek/10.0/toolbox/r2017aom'          % Ubuntu
% addpath 'C:\Program Files\Mosek\10.0\toolbox\r2017aom'    % Windows

close all
clear
rng('default')
warning('off','all')

%% Specify system model and simulation
% System model
system = 'quad2D';
define_system

%% Train GP models
% Specify the number of random datasets to draw and average over.
n_trials = 10;

Ts_max_vec = zeros(length(N_vec),n_trials);

%% Train models and compute minimum control frequency
for p = 1:n_trials
    fprintf("Trial number %.f \n", p);
    % Create noise-free training samples
    switch system
        case 'quad1D'
            [z_obs,y_obs] = quad1D_samples(N_vec(end),params);
        case 'quad2D'
            [z_obs,y_obs] = quad2D_samples(N_vec(end),params);
    end

    % Add observation noise
    for i = 1:n
        y_obs(:,i) = y_obs(:,i) + sigma_n(i)*randn(N_vec(end),1);
    end

    % Gradually increase dataset and simulate each time
    for j = 1:length(N_vec)
        % Fit GP models
        gpr = cell(n,1);
        for i = 1:n
            % Take subset of training samples (gradually increasing dataset)
            gpr{i} = fitrgp(z_obs(1:N_vec(j),:),y_obs(1:N_vec(j),i),'KernelFunction',...
                'ardsquaredexponential','Sigma',sigma_n(i),'ConstantSigma',true);
        end

        % Initialize derivative predictions
        dmu = zeros(n+m,n);
        dSigma = zeros(n+m,n+m,n);

        % Compute derivatives
        for i = 1:n
            [sigma_f,L] = get_hyperparameters(gpr{i});
            L_inv = L^(-2);
            Z_train = gpr{i}.X;
            K_matrix = get_K(Z_train,sigma_f,L_inv,sigma_n(i));
            dmu(:,i) = get_dmu(sigma_f,L_inv,K_matrix,y_obs(1:N_vec(j),i),Z_train,[x_s,u_s]);
            dSigma(:,:,i) = get_dSigma(sigma_f,L_inv,K_matrix,Z_train,[x_s,u_s]);
        end

        [A,B] = dmu2AB(dmu,n);
        [Au,Bu] = dSigma2AB(dSigma,n);

        % Transform the matrix-valued uncertainty
        [H,E,F] = transform_uncertainty(Au,Bu);

        % Maximize the sampling time, i.e., minimize the control frequency.
        eps_vec = logspace(-2,2,10);
        [Ts_max, K, Ts_vec] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec);

        % Store the maximum sampling interval
        Ts_max_vec(j,p) = Ts_max;
    end
end

% Plot minimum control frequency over the number of datapoints N
success_rate = zeros(1,length(N_vec));
if length(N_vec) > 1
    % Remove results where Ts_max = 0
    fc_min_mean = zeros(1,length(N_vec));
    fc_min_std = zeros(1,length(N_vec));
    for k = 1:length(N_vec)
        tmp1 = Ts_max_vec(k,:);
        tmp2 = tmp1(tmp1 >= 0.001);
        fc_min_mean(k) = mean(1./tmp2);
        fc_min_std(k) = std(1./tmp2);
        success_rate(k) = length(tmp2) / length(tmp1);
    end    
  
    % Plot minimum control frequency
    figure()
    errorbar(N_vec,fc_min_mean,fc_min_std);
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
