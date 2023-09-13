% Simulate the system several times for the same initial condition, but
% with different control frequencies and controllers.

% addpath '/home/ralf/mosek/10.0/toolbox/r2017aom'          % Ubuntu
% addpath 'C:\Program Files\Mosek\10.0\toolbox\r2017aom'    % Windows

close all
clear
warning('off','all')
colors = matlab_plot_colors;

% Where to load the data from
load_path = "./eval/data/analyze_performance/main_figure_v1/";

%% Specify system simulation
% System model
system = 'quad2D';
define_system

% Simulation parameters
Tsim = 10;
steps_per_sample = 11;
tspan_plot = linspace(0,Tsim,101);

%% Load learned system models and controllers
load(load_path + "A_all",'A_all');
load(load_path + "B_all",'B_all');
load(load_path + "Au_all",'Au_all');
load(load_path + "Bu_all",'Bu_all');
load(load_path + "N_vec",'N_vec');
load(load_path + "K_all",'K_all');
load(load_path + "fc_vec",'fc_vec');

n_trials = size(A_all,4);

%% Simulate systems
J_vec = zeros(length(n_trials), length(N_vec), length(fc_vec));
Ts_max_vec_only_feasible = cell(1,length(N_vec));
% zeros(length(N_vec),size(delta_x0,1),n_trials,length(xi));
for j = 1:length(N_vec)
    for q = 1:length(fc_vec)
        for p = 1:n_trials
            % Get system matrices and uncertainty
            A = A_all(:,:,j,p);
            B = B_all(:,:,j,p);
            Au = Au_all(:,:,j,p);
            Bu = Bu_all(:,:,j,p);
            K = K_all(:,:,j,p,q);

            % Check if there is a stabilizing controller for the given
            % model and frequency
            if all(K == 0)
                J_vec(p,j,q) = -1;
                continue;
            end

            % Transform uncertainty representation
            [H,E,E_u] = transform_uncertainty(Au,Bu);

            % Set sampling time
            Ts = 1/fc_vec(q);

            %% Simulate System
            J_tmp = zeros(1,size(delta_x0,1));
            % Use different initial states
            for r = 1:size(delta_x0,1)
                switch system
                    case 'quad1D'
                        x_k = x_s + delta_x0(r,:);
                    case 'quad2D'
                        x_k = x_s + delta_x0(r,:);
                end

                % Simulate system
                x_sim = []; u_sim = [];
                tspan = linspace(0,Ts*ceil(Tsim/Ts),(steps_per_sample-1)*ceil(Tsim/Ts)+1);
                for k = 1:ceil(Tsim/Ts)
                    % Calculate input, which is kept constant
                    u = K * (x_k - x_s)' + u_s';

                    [t,x_new] = ode45(@(t,x) f(x,u,params), linspace(0,Ts,steps_per_sample), x_k);
                    if k > 1
                        x_new = x_new(2:end,:);
                    end
                    x_sim = [x_sim; x_new];
                    u_sim = [u_sim; repmat(u',size(x_new,1),1)];

                    % Update x_k
                    x_k = x_new(end,:);
                end

                J_tmp(r) = compute_cost(x_sim,u_sim,x_s,u_s,J1,J2,Tsim);
            end
            % Calculate cost
            J_vec(p,j,q) = mean(J_tmp);
        end
    end
end

%% Scatter costs over N and fc
N_plot_feasible = [];
fc_plot_feasible = [];
J_plot = [];
N_plot_infeasible = [];
fc_plot_infeasible = [];

J_mean_mat = zeros(length(N_vec),length(fc_vec));
for j=1:length(N_vec)
    for q = 1:length(fc_vec)
        % Check if more than 50% of the trials were infeasible
        tmp1 = J_vec(:,j,q);
        tmp2 = tmp1(tmp1 ~= -1);
        % if length(tmp2) < 0.5 * length(tmp1)
        if length(tmp2) <= 2
            J_mean_mat(j,q) = -1;
            N_plot_infeasible = [N_plot_infeasible; N_vec(j)];
            fc_plot_infeasible = [fc_plot_infeasible; fc_vec(q)];
        else
            J_mean_mat(j,q) = mean(tmp2);
            N_plot_feasible = [N_plot_feasible; N_vec(j)];
            fc_plot_feasible = [fc_plot_feasible; fc_vec(q)];
            J_plot = [J_plot; mean(tmp2)];
        end
    end
end

% Draw colormap corresponding to cost values
scatter(N_plot_feasible,fc_plot_feasible,200,J_plot,'filled','Marker','square'); hold on
scatter(N_plot_infeasible,fc_plot_infeasible,200,'filled','white','Marker','square');
J_mean_mat(J_mean_mat == -1) = NaN;
% imagesc([N_vec(1), N_vec(end)], [fc_vec(1), fc_vec(end)], J_mean_mat,'AlphaData',~isnan(J_mean_mat)); hold on
colormap turbo
colorbar
% set(gca,'YDir','normal')

% Draw contour lines
[N_plot,fc_plot] = meshgrid(N_vec,fc_vec);
contour(N_plot,fc_plot,J_mean_mat',[5:10]);
% contour(N_plot,fc_plot,J_mean_mat',[5.5, 7, 8.5, 10]);
% contour(N_plot,fc_plot,J_mean_mat',[5.5, 7.5, 9.5]);

xlim([150 1050]);
ylim([9 31]);
cleanfigure();
