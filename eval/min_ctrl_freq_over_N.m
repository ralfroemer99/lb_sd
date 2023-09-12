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
n_trials = 4;
Ts_max_vec = zeros(length(N_vec),n_trials);

%% Train models and compute minimum control frequency
disp('Starting model learning')

for p = 1:n_trials
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

        % Save system matrices
        A_all(:,:,j,p) = A;         B_all(:,:,j,p) = B;
        Au_all(:,:,j,p) = Au;       Bu_all(:,:,j,p) = Bu;

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
if length(N_vec) > 1
    % Remove results where Ts_max = 0
    fc_min_mean = zeros(1,length(N_vec));
    fc_min_std = zeros(1,length(N_vec));
    for k = 1:length(N_vec)
        tmp = Ts_max_vec(k,:);
        tmp = tmp(tmp >= 0.01);
        fc_min_mean(k) = mean(1./tmp);
        fc_min_std(k) = std(1./tmp);
    end    
    figure()
    errorbar(N_vec,fc_min_mean,fc_min_std);
    xlabel('$N$','Interpreter','latex');
    ylabel('$f_{c,min}$','Interpreter','latex');
end