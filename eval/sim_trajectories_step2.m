% Simulate the system several times for the same initial condition, but
% with different control frequencies and controllers.

% addpath '/home/ralf/mosek/10.0/toolbox/r2017aom'          % Ubuntu
% addpath 'C:\Program Files\Mosek\10.0\toolbox\r2017aom'    % Windows

close all
clear
rng('default')
warning('off','all')
colors = matlab_plot_colors;

% Where to load the data from
load_path = "./eval/data/sim_trajectories/gam4p5/";

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
load(load_path + "Ts_max_vec",'Ts_max_vec');
load(load_path + "K_all",'K_all');
load(load_path + "xi",'xi');
n_trials = size(A_all,4);

%% Simulate systems
x_sim_plot = zeros(length(tspan_plot),length(N_vec),n_trials,length(xi),size(delta_x0,1));
z_sim_plot = zeros(length(tspan_plot),length(N_vec),n_trials,length(xi),size(delta_x0,1));
J_vec = cell(1,length(N_vec));
Ts_max_vec_only_feasible = cell(1,length(N_vec));
for j = 1:length(N_vec)
    J_tmp = zeros(size(delta_x0,1),sum(Ts_max_vec(j,:)~=0),length(xi));
    for q = 1:length(xi)
        tmp1 = 1;
        for p = 1:n_trials
            % Get system matrices and uncertainty
            A = A_all(:,:,j,p);
            B = B_all(:,:,j,p);
            Au = Au_all(:,:,j,p);
            Bu = Bu_all(:,:,j,p);
            K = K_all(:,:,j,p,q);

            % Transform uncertainty representation
            [H,E,E_u] = transform_uncertainty(Au,Bu);

            % Set sampling time
            Ts_max = Ts_max_vec(j,p);
            if Ts_max == 0
                continue
            end
            Ts = xi(q) * Ts_max;
            Ts_max_vec_only_feasible{j}(tmp1,q) = Ts;

            %% Simulate System
            % Initial state
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

                % Calculate cost
                J_tmp(r,tmp1,q) = compute_cost(x_sim,u_sim,x_s,u_s,J1,J2,Tsim);
                % J_vec(j,r,p,q) = compute_cost(x_sim,u_sim,x_s,u_s,J1,J2,Tsim);

                % Store trajectories
                x_sim_plot(:,j,p,q,r) = interp1(tspan,x_sim(:,1),tspan_plot)';
                z_sim_plot(:,j,p,q,r) = interp1(tspan,x_sim(:,3),tspan_plot)';
            end
            tmp1 = tmp1 + 1;
        end
    end
    J_vec{j} = J_tmp;
end

%% Plot
% Mean and standard deviation or min/max trajectories
plot_with_std = 1;

% Plot trajectories
which_x0 = 1;       % For which initial state
figure();
for q = 1:length(xi)
    legend_entries = cell(1,length(N_vec));
    for j = 1:length(N_vec)
        % Only consider runs where Ts_max > 0
        which_indices = find(Ts_max_vec(j,:) > 0);

        % Calculate mean trajectory
        x_sim_mean = mean(x_sim_plot(:,j,which_indices,q,which_x0),3)';
        z_sim_mean = mean(z_sim_plot(:,j,which_indices,q,which_x0),3)';

        % Get area to fill
        if plot_with_std
            x_sim_std = std(x_sim_plot(:,j,which_indices,q,which_x0),[],3)';
            [tmp1, tmp2] = shaded_plot_mean_std(tspan_plot, x_sim_mean, x_sim_std);
            z_sim_std = std(z_sim_plot(:,j,which_indices,q,which_x0),[],3)';
            [tmp3, tmp4] = shaded_plot_mean_std(tspan_plot, z_sim_mean, z_sim_std);
        else
            x_sim_min = min(x_sim_plot(:,j,which_indices,q,which_x0),[],3);
            x_sim_max = max(x_sim_plot(:,j,which_indices,q,which_x0),[],3);
            [tmp1, tmp2] = shaded_plot_min_max(tspan_plot, x_sim_min, x_sim_max);
            z_sim_min = min(z_sim_plot(:,j,which_indices,q,which_x0),[],3);
            z_sim_max = max(z_sim_plot(:,j,which_indices,q,which_x0),[],3);
            [tmp3, tmp4] = shaded_plot_min_max(tspan_plot, z_sim_min, z_sim_max);
        end

        % Plot x
        subplot(2,length(xi),2*q-1); ylabel('$x(t)$','interpreter','latex');
        fill(tmp1, tmp2, colors(j,:),'EdgeColor',colors(j,:),'FaceAlpha',0.2,'EdgeAlpha',0.2,'HandleVisibility','off'); hold on
        plot(tspan_plot,x_sim_mean, 'Color', colors(j,:), 'LineWidth', 1);

        % Plot z
        subplot(2,length(xi),2*q); ylabel('$z(t)$','interpreter','latex'); xlabel('$t$','interpreter','latex');
        fill(tmp3,tmp4,colors(j,:),'EdgeColor',colors(j,:),'FaceAlpha',0.2,'EdgeAlpha',0.2,'HandleVisibility','off'); hold on
        plot(tspan_plot,z_sim_mean,'Color',colors(j,:),'LineWidth', 1);

        % Legend
        tmp = Ts_max_vec(j,:);
        legend_entries{j} = "$N = " + num2str(N_vec(j)) + "$, $f_\mathrm{c} = " + num2str(...
            round(1/xi(q) * mean(1./tmp(tmp ~= 0)),1)) + "$";
    end
    legend(legend_entries, 'Interpreter','latex')
    sgtitle("$f_\mathrm{c} = " + num2str(1/xi(q)) + "f_\mathrm{c,min}$",'Interpreter','latex')
end

figure();
for j = 1:length(N_vec)
    legend_entries = cell(1,length(N_vec));
    for q = 1:length(xi)
        % Only consider runs where Ts_max > 0
        which_indices = find(Ts_max_vec(j,:) > 0);

        % Calculate mean trajectory
        x_sim_mean = mean(x_sim_plot(:,j,which_indices,q,which_x0),3)';
        z_sim_mean = mean(z_sim_plot(:,j,which_indices,q,which_x0),3)';

        % Get area to fill
        if plot_with_std
            x_sim_std = std(x_sim_plot(:,j,which_indices,q,which_x0),[],3)';
            [tmp1, tmp2] = shaded_plot_mean_std(tspan_plot, x_sim_mean, x_sim_std);
            z_sim_std = std(z_sim_plot(:,j,which_indices,q,which_x0),[],3)';
            [tmp3, tmp4] = shaded_plot_mean_std(tspan_plot, z_sim_mean, z_sim_std);
        else
            x_sim_min = min(x_sim_plot(:,j,which_indices,q,which_x0),[],3);
            x_sim_max = max(x_sim_plot(:,j,which_indices,q,which_x0),[],3);
            [tmp1, tmp2] = shaded_plot_min_max(tspan_plot, x_sim_min', x_sim_max');
            z_sim_min = min(z_sim_plot(:,j,which_indices,q,which_x0),[],3);
            z_sim_max = max(z_sim_plot(:,j,which_indices,q,which_x0),[],3);
            [tmp3, tmp4] = shaded_plot_min_max(tspan_plot, z_sim_min', z_sim_max');
        end

        % Plot x with standard deviation
        subplot(2,length(N_vec),j);
        if j == 1
            ylabel('$x(t)$','interpreter','latex');
        end
        fill(tmp1, tmp2, colors(q,:),'EdgeColor',colors(q,:),'FaceAlpha',0.2,'EdgeAlpha',0.2,'HandleVisibility','off'); hold on
        plot(tspan_plot,x_sim_mean, 'Color', colors(q,:), 'LineWidth', 1);
        title("$N = " + num2str(N_vec(j)) + "$",'Interpreter','latex')

        % Plot z with standard deviation
        subplot(2,length(N_vec),j+length(N_vec));
        if j == 1
            ylabel('$z(t)$','interpreter','latex');
        end
        xlabel('$t$','interpreter','latex');
        fill(tmp3,tmp4,colors(q,:),'EdgeColor',colors(q,:),'FaceAlpha',0.2,'EdgeAlpha',0.2,'HandleVisibility','off'); hold on
        plot(tspan_plot,z_sim_mean,'Color',colors(q,:),'LineWidth', 1);

        tmp = Ts_max_vec(j,:);
        legend_entries{q} = "$f_\mathrm{c} = " + num2str(round(1/xi(q) * mean(1./tmp(tmp ~= 0)),1)) + "$";
    end
    legend(legend_entries, 'Interpreter', 'latex')
    % sgtitle("$N = " + num2str(N_vec(i)) + "$",'Interpreter','latex')
end

% Box plot of the costs
figure()
x_pos = [0.25 0.5 0.75;
         1.25 1.5 1.75;
         2.25 2.5 2.75];
for i = 1:length(N_vec)
    J_tmp = [];
    % labels = [];
    for j = 2:length(xi)
        J_tmp = [J_tmp, reshape(J_vec{i}(:,:,j),[],1)];
        % labels = [labels; repmat({"N = " + num2str(N_vec(i))},length(reshape(J_vec{i}(:,:,j),[],1)),1)];
    end
    boxplot(J_tmp,'positions',x_pos(:,i)','Widths',0.2,'Colors',colors(i,:)); hold on
end
