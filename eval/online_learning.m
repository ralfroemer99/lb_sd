% addpath '/home/ralf/mosek/10.0/toolbox/r2017aom'          % Ubuntu
% addpath 'C:\Program Files\Mosek\10.0\toolbox\r2017aom'    % Windows
close all
clear
rng('default')

%% Define system model
% Real system
% A = [0 1 0; 0 -3 10; 0 0 -5];
% B = [0; 0; 1];
x_s = [1,0,0,0,0,0]; u_s = [0.981,0.981];
[A,B] = quad2D_linearization(x_s,u_s,quad2D_params);
[n,m] = size(B);

% Noise variance
sigma_n = 1;

% Specify after how many data points to repeat the offline phase
N_total = 400;
% dN = 100;
N_offline = [50,150,N_total];
f_smaller = 1:-0.1:0.4;

% State and input bounds for generating data
xlim_low = -5*ones(n,1);    xlim_high = 5*ones(n,1);
ulim_low = -5*ones(m,1);    ulim_high = 5*ones(m,1);

% Draw random training data
N_init = 80;
[x_rand,u_rand] = state_input_samples(xlim_low,xlim_high,ulim_low,ulim_high,N_total);

%% Initialization
% Define initial estimate of the model
A0 = zeros(n);
SigmaA0 = cell(1,n);
for i = 1:n
    SigmaA0{i} = 5*eye(n);
end

% Store stuff
alpha_Ai_vec = [zeros(1,N_total)];         % Uncertainty magnitude of i-th estimate
alpha_Ei_vec = [zeros(1,N_total)];         % Uncertainty magnitude of interval centered at A
alpha_vec = zeros(length(f_smaller),length(N_offline));
fc0_vec = zeros(1,length(N_offline));

i = 1;
for iter = 1:length(N_offline)
    %% Offline phase
    if iter == 1
        % Observe initial data
        for i = 1:N_init
            x = x_rand(:,i);
            u = u_rand(:,i);
            dx = A * x + B * u + sigma_n * rand(n,1);
            [A0, SigmaA0] = blr_update_A(x,u,dx,A0,B,SigmaA0,sigma_n);
        end
    end

    % Get confidence interval for A0
    DeltaA0 = Sigma2confidenceInterval(SigmaA0);

    % Compute initial control frequency
    Ts0 = max_Ts_additive_uncertainty(A0,B,sqrt(n)*eye(n),induced_norm12(DeltaA0));
    if Ts0 == 0
        error('Need more initial training data!')
    end
    fc0_vec(iter) = 1/Ts0;

    alpha_vec(1,iter) = induced_norm12(DeltaA0);

    % Compute maximum uncertainty bound for smaller frequencies
    % F = f_smaller*fc0_vec(iter);
    for k = 2:length(f_smaller)
        [alpha_vec(k,iter),K_tmp] = max_disturbance(A0,B,sqrt(n)*eye(n),1/(f_smaller(k)*fc0_vec(iter)));
    end

    %% Online phase
    % Aiold = A0;
    % SigmaAiold = SigmaA0;
    % DeltaAiold = DeltaA0;
    Ai = A0;
    SigmaAi = SigmaA0;

    % iter_jump = 50;
    % iscontained_vec = [];
    if iter == 1
        N_start = 1;
    else
        N_start = N_offline(iter-1) + 1;
    end
    for i = N_start:N_offline(iter)
        % Collect data
        x = x_rand(:,i);
        u = u_rand(:,i);

        % Update model
        dx = A * x + B * u + sigma_n * rand(n,1);
        [Ai, SigmaAi] = blr_update_A(x,u,dx,Ai,B,SigmaAi,sigma_n);
        DeltaAi = Sigma2confidenceInterval(SigmaAi);

        % --> Compute interval centered at A0 that contains Ai +- DeltaAi
        Ci = A0 - (Ai - DeltaAi);
        Di = A0 - (Ai + DeltaAi);
        Ei = max(abs(Ci),abs(Di));

        % Store uncertainty magnitudes of Ai and Ei
        alpha_Ai_vec(i) = induced_norm12(DeltaAi);
        alpha_Ei_vec(i) = induced_norm12(Ei);

        % if mod(i,iter_jump) == 0
        %     % Check if the control frequency can be reduced
        % 
        %     % Check if new uncertainty interval is contained in previous one
        %     iscontained_vec(end+1) = all((Aiold - DeltaAiold) <= (Ai - DeltaAi) ...
        %         & (Aiold + DeltaAiold) >= (Ai - DeltaAi),'all');
        % 
        %     % Store estimates for next comparison
        %     Aiold = Ai;
        %     SigmaAiold = SigmaAi;
        %     DeltaAiold = DeltaAi;
        % end
    end

    % Update model
    A0 = Ai;
    SigmaA0 = SigmaAi;
end

plot(1:N_total,alpha_Ai_vec); hold on
plot(1:N_total,alpha_Ei_vec);
xlabel('$i$','Interpreter','latex');
ylabel('Uncertainty')
% legend('$\|\mathbf{A}_i\|_{1,2}$','$\|\mathbf{E}_i\|_{1,2}$','Interpreter','latex')

% Draw horizontal lines
% figure
% plot(0:iter_jump:N_iter,alpha_Ei_vec(1:iter_jump:end),'b'); hold on
% plot(0:iter_jump:N_iter,alpha_Ai_vec(1:iter_jump:end),'b--');
% plot([0 N_iter], [induced_norm12(DeltaA0),induced_norm12(DeltaA0)])
% legendentries = {'$\|\mathbf{\tilde{\Delta}}_i^A\|_{1,2}$','$\|\mathbf{\Delta}_i^A\|_{1,2}$','$\alpha^*(f_{\mathrm{c},0})$'};
legendentries = {'$\|\mathbf{\tilde{\Delta}}_i^A\|_{1,2}$','$\|\mathbf{\Delta}_i^A\|_{1,2}$'};
for i = 1:length(N_offline)
    if i == 1
        N_start = 1;
    else
        N_start = N_offline(i-1) + 1;
    end
    for k = 1:length(f_smaller)
        if alpha_vec(k,i) ~= 0
            plot([N_start N_offline(i)], [alpha_vec(k,i),alpha_vec(k,i)]);
            % legendentries{end+1} = "$\alpha^*(" + num2str(f_smaller(k)/fc0_vec(i)) + "f_{\mathrm{c},0})$";
            legendentries{end+1} = "$f_\mathrm{c} = " + num2str(f_smaller(k)*fc0_vec(i)) + "$";
        end
    end
end

xlabel('$i$','Interpreter','latex');
ylabel('Uncertainty')
legend(legendentries,'Interpreter','latex')
