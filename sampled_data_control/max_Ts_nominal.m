function [Ts_max,K] = max_Ts_nominal(A,B)
small_scalar = 1e-8;

% Get dimensions
[n,m] = size(B);

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

% Stability LMIs
lmi1_lhs = ...
    [Z1, Z2, Q2';
     Z2', Z3, Q3';
     Q2, Q3, -R];
lmi1_rhs = -Ts_inv * ...
    [Q2+Q2', Q3-Q2'+Q1*A'+Y'*B', zeros(n,n);
     Q3'-Q2+A*Q1'+B*Y, -Q3-Q3', zeros(n,n);
     zeros(n), zeros(n), zeros(n)];
lmi1 = lmi1_lhs <= lmi1_rhs - small_scalar*eye(3*n);
lmi2 = ...
    [2*Q1-R, zeros(n), Y'*B';
     zeros(n), Z1, Z2;
     B*Y, Z2, Z3] >= 0;
const = [Ts_inv >= 1, Q1 >= small_scalar*eye(n), R >= small_scalar*eye(n), lmi1,lmi2, 1000 >= Ts_inv >= 1];

% Objective
obj = Ts_inv;

% Solve using the bisection method
diagnostics = bisection(const, obj, sdpsettings('solver','mosek','verbose',0));

% Get sampling time and control gain
Ts_max = 1/value(Ts_inv);
K = value(Y)/value(Q1);

end
