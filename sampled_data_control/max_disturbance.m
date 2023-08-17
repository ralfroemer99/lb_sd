function [alpha_max,K] = max_disturbance(A,B,H,Ts_max)
eps = 1e-8;
Ts_inv = 1/Ts_max;

[n,m] = size(B);
% Define optimization variables

gam = sdpvar(1,1);
Q1 = sdpvar(n);           % symmetric
Q2 = sdpvar(n,n);
Q3 = sdpvar(n,n);
Z1 = sdpvar(n,n);
Z2 = sdpvar(n,n);
Z3 = sdpvar(n,n);
Rbar = sdpvar(n);         % symmetric
Y = sdpvar(m,n);

% Stability LMIs
lmi1_lhs = [Z1, Z2, zeros(n,n), Q2', zeros(n,n);
            Z2', Z3, zeros(n,n), Q3', zeros(n,n);
            zeros(n,n), zeros(n,n), zeros(n,n), zeros(n,n), zeros(n,n);
            Q2, Q3, zeros(n,n), -Rbar, zeros(n,n);
            zeros(n,n), zeros(n,n), zeros(n,n), zeros(n,n), zeros(n,n)];
lmi1_rhs = -Ts_inv*...
    [Q2+Q2', Q1*A'+Y'*B'-Q2'+Q3, zeros(n,n), zeros(n,n), Q1*H';
     A*Q1+B*Y-Q2+Q3', -Q3-Q3', eye(n), zeros(n,n), zeros(n,n);
     zeros(n,n), eye(n), -eye(n), zeros(n,n), zeros(n,n);
     zeros(n,n), zeros(n,n), zeros(n,n), zeros(n,n), zeros(n,n);
     H*Q1, zeros(n,n), zeros(n,n), zeros(n,n), -gam*eye(n)];
lmi1 = (lmi1_lhs <= lmi1_rhs - eps*eye(5*n));

lmi2 = [2*Q1-Rbar, zeros(n,n), Y'*B';
        zeros(n,n), Z1, Z2;
        B*Y, Z2', Z3] >= 0;

const = [Q1 >= eps*eye(n), Rbar >= eps*eye(n),lmi1,lmi2];

% Objective
obj = gam;

% Solve using the bisection method
diagnostics = bisection(const, obj, sdpsettings('solver','mosek','verbose',0));

% Get sampling time and control gain
if diagnostics.problem == 0
    alpha_max = 1 / sqrt(value(gam));
    K = value(Y)/value(Q1);
else
    alpha_max = 0;
    K = zeros(m,n);
end

end

