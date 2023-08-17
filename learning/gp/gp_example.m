close all
clear
%% Initialize and train GP models
rng('default') % For reproducibility
sigma_n = 1;
x_observed2 = 0 + 10*rand(30,1);
x_observed1 = x_observed2(1:15);
% x_observed1 = linspace(0,10,11)';
% x_observed2 = linspace(0,10,31)';
y_observed1 = x_observed1.*sin(x_observed1) + 0.2*randn(size(x_observed1));
y_observed2 = x_observed2.*sin(x_observed2) + 0.2*randn(size(x_observed2));

gpr1 = fitrgp(x_observed1,y_observed1,'Sigma',sigma_n,'ConstantSigma',true);
gpr2 = fitrgp(x_observed2,y_observed2,'Sigma',sigma_n,'ConstantSigma',true);

x = linspace(0,10,51)';
[mu_pred1,sd_pred1] = predict(gpr1,x);
[mu_pred2,sd_pred2] = predict(gpr2,x);

%% Get derivative prediction
% Get kernel hyperparameters
[sigma_f1,L1] = get_hyperparameters(gpr1);
[sigma_f2,L2] = get_hyperparameters(gpr2);
L_inv1 = L1^(-2);
L_inv2 = L2^(-2);
Z_train1 = gpr1.X;
Z_train2 = gpr2.X;
K1 = get_K(Z_train1,sigma_f1,L_inv1,sigma_n);
K2 = get_K(Z_train2,sigma_f2,L_inv2,sigma_n);
dmu_pred1 = zeros(length(x),1);
dmu_pred2 = zeros(length(x),1);
dSigma_pred1 = zeros(length(x),1);
dSigma_pred2 = zeros(length(x),1);
for i=1:length(x)
    dmu_pred1(i) = get_dmu(sigma_f1,L_inv1,K1,y_observed1,Z_train1,x(i));
    dSigma_pred1(i) = get_dSigma(sigma_f1,L_inv1,K1,Z_train1,x(i));
    dmu_pred2(i) = get_dmu(sigma_f2,L_inv2,K2,y_observed2,Z_train2,x(i));
    dSigma_pred2(i) = get_dSigma(sigma_f2,L_inv2,K2,Z_train2,x(i));
end

%% Plot
fig = figure;

tiledlayout(2,2,'TileSpacing','compact')

% Prediction
nexttile
hold on
scatter(x_observed1,y_observed1,'r') % Observed data points
plot(x,x.*sin(x),'--r')
plot(x,mu_pred1,'b','LineWidth',1)                  % GPR predictions
patch([x;flipud(x)],[mu_pred1 - 2*sd_pred1;flipud(mu_pred1 + 2*sd_pred1)],'k','FaceAlpha',0.1); % Prediction intervals
hold off
box on
title('GP prediction')
legend({'Training data','$g(z) = z\sin{(z)}$','$\mu(z)$'},'Interpreter','latex')

% Derivative prediction
nexttile
hold on
plot(x,x.*cos(x) + sin(x),'--r')         
plot(x,dmu_pred1,'b','LineWidth',1)
patch([x;flipud(x)],[dmu_pred1 - 2*sqrt(dSigma_pred1);flipud(dmu_pred1 + 2*sqrt(dSigma_pred1))],'k','FaceAlpha',0.1); % Prediction intervals
hold off
box on
title('GP derivative prediction')
legend({'$dg(z) = \sin{(z)}+z\cos{(z)}$','$d\mu(z)$'},'Interpreter','latex')

% Prediction
nexttile
hold on
scatter(x_observed2,y_observed2,'r') % Observed data points
plot(x,x.*sin(x),'--r')
plot(x,mu_pred2,'b','LineWidth',1)                   % GPR predictions
patch([x;flipud(x)],[mu_pred2 - 2*sd_pred2;flipud(mu_pred2 + 2*sd_pred2)],'k','FaceAlpha',0.1); % Prediction intervals
hold off
box on
% legend({'Training data','$g(z) = z\sin{(z)}$','$\mu(z)$'},'Interpreter','latex')

% Derivative prediction
nexttile
hold on
plot(x,x.*cos(x) + sin(x),'--r')
plot(x,dmu_pred2,'b','LineWidth',1)
patch([x;flipud(x)],[dmu_pred2 - 2*sqrt(dSigma_pred2);flipud(dmu_pred2 + 2*sqrt(dSigma_pred2))],'k','FaceAlpha',0.1); % Prediction intervals
hold off
box on
% legend({'$dg(z) = \sin{(z)}+z\cos{(z)}$','$d\mu(z)$'},'Interpreter','latex')
