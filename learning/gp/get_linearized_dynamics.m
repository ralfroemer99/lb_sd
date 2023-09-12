function [A,B,Au,Bu] = get_linearized_dynamics(gpr,sigma_n,y_obs,x_s,u_s)
% Compute the uncertain linearization of a learned GP dynamics model at the
% state-input pair (x_s,u_s). The linearization has the form
% dx = (A + Au o Omega) * x + (B + Bu o Psi) * u
% where Omega (Psi) is an n x n (n x m) matrix with unknown elements in 
% [-1,1] and o denotes the element-wise product.

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
