function [A,B,Au,Bu] = get_linearized_dynamics(gpr,sigma_n,y_obs,x_s,u_s)
% Get dimensions
n = length(x_s);
m = length(u_s);

% Initialize derivative predictions
dmu = zeros(n+m,n);
dSigma = zeros(n+m,n+m,n);

% Iterate
for i = 1:n
    [sigma_f,L] = get_hyperparameters(gpr{i});
    L_inv = L^(-2);
    Z_train = gpr{i}.X;
    K_matrix = get_K(Z_train,sigma_f,L_inv,sigma_n(i));
    dmu(:,i) = get_dmu(sigma_f,L_inv,K_matrix,y_obs(:,i),Z_train,[x_s,u_s]);
    dSigma(:,:,i) = get_dSigma(sigma_f,L_inv,K_matrix,Z_train,[x_s,u_s]);
end

% Linearized uncertain system
[A,B] = dmu2AB(dmu,n);
[Au,Bu] = dSigma2AB(dSigma,n);
end
