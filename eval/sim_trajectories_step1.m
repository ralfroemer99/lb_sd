% Learn the nonlinear dynamics from measured data and design an optimized
% sampled-data controller for the linearized system.

% addpath '/home/ralf/mosek/10.0/toolbox/r2017aom'          % Ubuntu
% addpath 'C:\Program Files\Mosek\10.0\toolbox\r2017aom'    % Windows

close all
clear
rng('default')
warning('off','all')

% Specify save path
save_path = "./eval/data/sim_trajectories/";

%% Specify system model and simulation
% System model
system = 'quad2D';
define_system

% Specify the factors between the used sampling time and the maximum
% sampling interval for which to compute controllers.
xi = 1./[1, 1.25, 1.5, 2];    % 1.25

%% Train GP models
n_trials = 5;
Ts_max_vec = zeros(length(N_vec),n_trials);

%% Train models and compute minimum control frequency
disp('Starting model learning')

% Storage variables for the system matrices
A_all = zeros(n,n,length(N_vec),n_trials);
B_all = zeros(n,m,length(N_vec),n_trials);
Au_all = zeros(n,n,length(N_vec),n_trials);
Bu_all = zeros(n,m,length(N_vec),n_trials);
K_all = zeros(m,n,length(N_vec),n_trials,length(xi));

% Compute the minimum control frequency.
for p = 1:n_trials
    % Create noise-free training samples
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

    % Gradually increase dataset and simulate each time
    for j = 1:length(N_vec)
        % Fit GP models
        gpr = cell(n,1);
        for i = 1:n
            % Take subset of training samples (growing dataset)
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

        % Transform uncertainty representation
        [H,E,F] = transform_uncertainty(Au,Bu);

        % Maximize sampling time
        eps_vec = logspace(-2,2,10);
        [Ts_max, K, Ts_vec] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec);

        % Store MSI
        Ts_max_vec(j,p) = Ts_max;

        % Save control gain
        K_all(:,:,j,p,1) = K;
    end
end

pause(1)
disp('Starting controller optimization')

%% Compute optimized controllers for different control frequencies
for q = 2:length(xi)
    fprintf("Optimizing controller for fc = %.2f * fc_min \n", 1/xi(q))
    for j = 1:length(N_vec)        % Don't optimize for all N
        for p = 1:n_trials
            % Get system matrices and uncertainty
            A = A_all(:,:,j,p);
            B = B_all(:,:,j,p);
            Au = Au_all(:,:,j,p);
            Bu = Bu_all(:,:,j,p);

            % Transform uncertainty representation
            [H,E,F] = transform_uncertainty(Au,Bu);

            % Set sampling time
            Ts_max = Ts_max_vec(j,p);
            if Ts_max == 0
                continue
            end
            Ts = xi(q) * Ts_max;

            % Maximize performance for a given sampling time
            eps_vec = logspace(-2,2,10);
            [eta, K, eta_vec] = min_J_norm_bounded(A,B,H,E,F,Ts,J1,J2,eps_vec);

            % Save control gain
            K_all(:,:,j,p,q) = K;
        end
    end
end

% Save data
save(save_path + "A_all",'A_all');
save(save_path + "B_all",'B_all');
save(save_path + "Au_all",'Au_all');
save(save_path + "Bu_all",'Bu_all');
save(save_path + "N_vec",'N_vec');
save(save_path + "Ts_max_vec",'Ts_max_vec');
save(save_path + "K_all",'K_all');
save(save_path + "xi",'xi');
