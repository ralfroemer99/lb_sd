function [Ts_max,K,Ts_vec] = max_Ts_norm_bounded(A,B,H,E,E_u,eps_vec)
small_scalar = 1e-8;

% Get dimensions
[n,m] = size(B);
p = size(E,1);

Ts_vec = zeros(length(eps_vec),1);
K_all = zeros(m,n,length(eps_vec));

% Stability LMIs --> Solve in parallel for different choices of eps
parfor i=1:length(eps_vec)
    eps2 = eps_vec(i);
    eps1 = 1/eps2;

    % Define optimization variables
    Ts_inv = sdpvar;
    Q1 = sdpvar(n,n);           % symmetric
    Q2 = sdpvar(n,n);
    Q3 = sdpvar(n,n);
    Z1 = sdpvar(n,n);
    Z2 = sdpvar(n,n);
    Z3 = sdpvar(n,n);
    R = sdpvar(n,n);         % symmetric
    Y = sdpvar(m,n);

    lmi1_lhs = ...
        [Z1, Z2, Q2', zeros(n,p), zeros(n,p);
         Z2', Z3, Q3', zeros(n,p), zeros(n,p);
         Q2, Q3, -R, zeros(n,p), zeros(n,p);
         zeros(p,n), zeros(p,n), zeros(p,n), zeros(p), zeros(p);
         zeros(p,n), zeros(p,n), zeros(p,n), zeros(p), zeros(p)];
    
    lmi1_rhs = -Ts_inv * ...
        [Q2+Q2', Q3-Q2'+Q1*A'+Y'*B', zeros(n,n), zeros(n,p), eps1*(Q1*E'+Y'*E_u');
         Q3'-Q2+A*Q1'+B*Y, -Q3-Q3', zeros(n,n), H, zeros(n,p);
         zeros(n), zeros(n), zeros(n), zeros(n,p), zeros(n,p);
         zeros(p,n), H', zeros(p,n), -eps1*eye(p), zeros(p);
         eps1*(E*Q1'+E_u*Y), zeros(p,n), zeros(p,n), zeros(p), -eps1*eye(p)];
    lmi1 = lmi1_lhs <= lmi1_rhs - small_scalar*eye(3*n + 2*p);  % Always infeasible
    
    lmi2 = ...
        [2*Q1-R, zeros(n), Y'*B', zeros(n,p), eps2*Y'*E_u';
         zeros(n), Z1, Z2, zeros(n,p), zeros(n,p);
         B*Y, Z2, Z3, H, zeros(n,p);
         zeros(p,n), zeros(p,n), H', eps2*eye(p), zeros(p);
         eps2*E_u*Y, zeros(p,n), zeros(p,n), zeros(p), eps2*eye(p)] >= 0;
    const = [100 >= Ts_inv >= 0.1, Q1 >= small_scalar*eye(n), R >= small_scalar*eye(n), lmi1, lmi2];

    % Objective
    obj = Ts_inv;
    
    % Solve using the bisection method
    diagnostics = bisection(const, obj, sdpsettings('solver','mosek','verbose',0));
    % disp(diagnostics.problem)

    % If feasible, check if solution is the best obtained so far
    
    if diagnostics.problem == 0
        Ts_vec(i) = 1/value(Ts_inv);
        K_all(:,:,i) = value(Y)/value(Q1);
    end
end

% Find best solution
Ts_max = max(Ts_vec);
idx = find(Ts_vec == Ts_max);
if length(idx) > 1
    K = K_all(:,:,idx(1));
else
    K = K_all(:,:,idx);
end

end
