% Learn the nonlinear dynamics from measured data. The obtained approximate
% dynamics are linearized and the norm-bounded uncertainty in A and B is
% computed. A sampled-data controller is designed for the linearized system
% such that the sampling time is maximal while ensuring stability of the
% system.

% addpath '/home/ralf/mosek/10.0/toolbox/r2017aom'          % Ubuntu
% addpath 'C:\Program Files\Mosek\10.0\toolbox\r2017aom'    % Windows
close all
clear
rng('default')

%% Specify system model
system = 'quad2D';
switch system
    case 'quad1D'
        n = 3;
        m = 1;
        z_s = [1,0,0,0];
        sigma_n = [0.5,0.5,1];
        N_vec = 60:20:200;
    case 'quad2D'
        n = 6;
        m = 2;
        z_s = [1,0,0,0,0,0,0.981,0.981];
        sigma_n = [0.1,0.1,0.1,0.1,0.1,0.1];
        % N_vec = 200:200:1000;
        N_vec = 1000;
        params = quad2D_params;
    otherwise
        error('Wrong system type!')
end

%% Initialize and train GP models
n_trials = 1;
Ts_max_vec = zeros(length(N_vec),n_trials);

for k=1:length(N_vec)
    for j = 1:n_trials
        N = N_vec(k);
        % Create noise-free samples
        switch system
            case 'quad1D'
                [z_obs,y_obs] = quad1D_samples(N);
            case 'quad2D'
                [z_obs,y_obs] = quad2D_samples(N,params);
        end
        
        % Add noise
        for i = 1:n
            y_obs(:,i) = y_obs(:,i) + sigma_n(i)*randn(N,1);
        end
        
        % Fit GP models
        gpr = cell(n,1);
        for i = 1:n
            gpr{i} = fitrgp(z_obs,y_obs(:,i),'KernelFunction','ardsquaredexponential',...
                'Sigma',sigma_n(i),'ConstantSigma',true);
        end
        
        %% Get derivative prediction
        % Initialize
        dmu = zeros(n+m,n);
        dSigma = zeros(n+m,n+m,n);
        
        % Compute derivatives
        for i = 1:n
            [sigma_f,L] = get_hyperparameters(gpr{i});
            L_inv = L^(-2);
            Z_train = gpr{i}.X;
            K = get_K(Z_train,sigma_f,L_inv,sigma_n(i));
            dmu(:,i) = get_dmu(sigma_f,L_inv,K,y_obs(:,i),Z_train,z_s);
            dSigma(:,:,i) = get_dSigma(sigma_f,L_inv,K,Z_train,z_s);
        end
        
        [A,B] = dmu2AB(dmu,n);
        [Au,Bu] = dSigma2AB(dSigma,n);
        disp('A = '); disp(A);
        disp('Au = '); disp(Au);
        disp('B = '); disp(B);
        disp('Bu = '); disp(Bu);
        
        %% Sampled-Data
        % Transform uncertainty representation
        [H,E,E_u] = transform_uncertainty(Au,Bu);
        
        % Maximize sampling time
        eps_vec = logspace(-3,3,10);
        [Ts_max, K, Ts_vec] = max_Ts_norm_bounded(A,B,H,E,E_u,eps_vec);
        fprintf('The maximum sampling time with uncertainty is: %.3f s. \n', Ts_max);
        [Ts_max_nom, K_nom] = max_Ts_nominal(A,B);
        fprintf('The maximum sampling time without uncertainty is: %.3f s. \n', Ts_max_nom);
    
        Ts_max_vec(k,j) = Ts_max;
    end
end

%% Plot
if length(N_vec) == 1
    % h = gca;
    % [X,Y] = meshgrid(eps_vec, eps_vec);
    % surf(X,Y,Ts_vec)
    % set(h,'Xscale','log'); set(h,'YScale','log')
    % xlabel('$\epsilon_3$','Interpreter','latex');
    % ylabel('$\epsilon_4$','Interpreter','latex');
    semilogx(eps_vec,Ts_vec);
    xlabel('$\epsilon_3$','Interpreter','latex');
end

if length(N_vec) > 1
    figure()
    errorbar(N_vec,mean(Ts_max_vec,2),2*std(Ts_max_vec,1,2));
    xlabel('$N$','Interpreter','latex');
    ylabel('$T_{s,max}$','Interpreter','latex');
    
    figure()
    f_min_vec = 1./Ts_max_vec;
    errorbar(N_vec,mean(f_min_vec,2),2*std(f_min_vec,1,2));
    xlabel('$N$','Interpreter','latex');
    ylabel('$f_{c,min}$','Interpreter','latex');
end

% matlab2tikz('plots/quad1D_fcmin.tex')

