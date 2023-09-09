function [eta_min,K,eta_vec] = min_J_norm_bounded(A,B,H,E,E_u,Ts,J1,J2,eps_vec)
small_scalar = 1e-8;

% Get dimensions
[n,m] = size(B);
p = size(E,1);

% Get all possible pairs of eps3 and eps4
[tmp1, tmp2] = meshgrid(eps_vec,eps_vec);
eps_mat = [tmp1(:) tmp2(:)];

% Create variables to store eta and K
eta_vec = zeros(1,length(eps_vec)^2);
K_all = zeros(m,n,length(eps_vec)^2);

% Stability LMIs
tmp = size(eps_mat,1);
parfor i = 1:tmp
    eps3 = eps_mat(i,1);
    eps4 = eps_mat(i,2);

    % Define optimization variables
    eta = sdpvar;
    Q1 = sdpvar(n,n);           % symmetric
    Q2 = sdpvar(n,n);
    Q3 = sdpvar(n,n);
    Z1 = sdpvar(n,n);
    Z2 = sdpvar(n,n);
    Z3 = sdpvar(n,n);
    R = sdpvar(n,n);         % symmetric
    Y = sdpvar(m,n);

    % Build matrix inequalities
    lmi1_lhs = ...
        [Z1, Z2, Q2', zeros(n,p), zeros(n,p);
        Z2', Z3, Q3', zeros(n,p), zeros(n,p);
        Q2, Q3, -R, zeros(n,p), zeros(n,p);
        zeros(p,n), zeros(p,n), zeros(p,n), zeros(p), zeros(p);
        zeros(p,n), zeros(p,n), zeros(p,n), zeros(p), zeros(p)];
    lmi1_rhs = -1/Ts * ...
        [Q2+Q2', Q3-Q2'+Q1*A'+Y'*B', zeros(n,n), zeros(n,p), eps3*(Q1*E'+Y'*E_u');
        Q3'-Q2+A*Q1'+B*Y, -Q3-Q3', zeros(n,n), H, zeros(n,p);
        zeros(n), zeros(n), zeros(n), zeros(n,p), zeros(n,p);
        zeros(p,n), H', zeros(p,n), -eps3*eye(p), zeros(p);
        eps3*(E*Q1'+E_u*Y), zeros(p,n), zeros(p,n), zeros(p), -eps3*eye(p)];
    lmi1 = lmi1_lhs <= lmi1_rhs - small_scalar*eye(3*n + 2*p);  % Always infeasible
    lmi2 = ...
        [2*Q1-R, zeros(n), Y'*B', zeros(n,p), eps4*Y'*E_u';
        zeros(n), Z1, Z2, zeros(n,p), zeros(n,p);
        B*Y, Z2, Z3, H, zeros(n,p);
        zeros(p,n), zeros(p,n), H', eps4*eye(p), zeros(p);
        eps4*E_u*Y, zeros(p,n), zeros(p,n), zeros(p), eps4*eye(p)] >= 0;
    lmi3_lhs = ...
        [zeros(n), Q1, Y';
        Q1', -inv(J1), zeros(n,m);
        Y, zeros(m,n), -inv(J2)];
    lmi3_rhs = eta * ...
        [eye(n), zeros(n), zeros(n,m);
        zeros(n), zeros(n), zeros(n,m);
        zeros(m,n), zeros(m,n), zeros(m)];
    lmi3 = lmi3_lhs <= lmi3_rhs;
    const = [eta >= small_scalar, Q1 >= small_scalar*eye(n), ...
        R >= small_scalar*eye(n), lmi1, lmi2, lmi3];

    % Objective
    obj = eta;

    % Solve using the bisection method
    diagnostics = bisection(const, obj, sdpsettings('solver','mosek','verbose',0));
    % disp(diagnostics.problem)

    % If feasible, check if solution is the best obtained so far
    if diagnostics.problem == 0
        eta_vec(i) = value(eta);
        K_all(:,:,i) = value(Y)/value(Q1);
    end
end

% Check if the problem was feasible
if all(eta_vec == 0)        % Infeasible: Return -1
    eta_min = -1;
    K = zeros(m,n);
else                        % Feasible: Return optimal controller
    eta_min = min(eta_vec(eta_vec > 0));
    idx = find(eta_vec == eta_min);
    K = K_all(:,:,idx);
end
end
