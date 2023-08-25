close all
clear

% Define state trajectory
T0 = 0;
tspan0 = 0:0.1:1;
tspan1 = 1:0.1:1.6;
tspan2 = 1.6:0.1:3;

% Define multiplication factors
k = [1, 1.5, 0.5];
tau = [2, 1.5, 4];

% Compute V_i(t)
V0 = k(1)*exp(-tspan0/tau(1));
V1 = k(2)*exp(-tspan1/tau(2));
V2 = k(3)*exp(-tspan2/tau(3));

% Compute V_i'(t)
V0_hat = V0;
gamma1 = V0_hat(end) / V1(1);
V1_hat = gamma1 * V1;
gamma2 = V1_hat(end) / V2(1);
V2_hat = gamma2 * V2;

% Plot V_i(t) and V_i'(t)
plot(tspan0, V0, 'b'); hold on
plot(tspan1, V1); 
plot(tspan2, V2); 
plot(tspan1, V1_hat, 'b'); 
plot(tspan2, V2_hat, 'b'); 


ylim([0 1])

matlab2tikz()
