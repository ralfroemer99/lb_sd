function is_feasible_mat = eval_controller_norm_bounded(A,B,H,E,E_u,K,Ts,eps_vec)
small_scalar = 1e-6;

A1 = B*K;
E1 = E_u*K;

[n,m] = size(B);
p = size(E,1);
% Define optimization variables
P1 = sdpvar(n,n);           % symmetric
P2 = sdpvar(n,n,'full');
P3 = sdpvar(n,n,'full');
Z1 = sdpvar(n,n,'full');
Z2 = sdpvar(n,n,'full');
Z3 = sdpvar(n,n,'full');
R = sdpvar(n,n);            % symmetric

eps1_vec = eps_vec;
eps2_vec = eps_vec;

is_feasible_mat = zeros(length(eps1_vec), length(eps2_vec));

% Stability LMIs
for i=1:length(eps1_vec)
    eps1 = eps1_vec(i);
    for j = 1:length(eps2_vec)
        eps2 = eps2_vec(j);
        P = [P1 zeros(n); P2 P3];
        Z = [Z1 Z2; Z2' Z3];
        W = P' * [zeros(n), eye(n); A+A1, -eye(n)] + ...
            [zeros(n), eye(n); A+A1, -eye(n)]' * P + ...
            Ts*Z + Ts*[zeros(n), zeros(n); zeros(n), Ts*R];
        lmi1 = [W, [P2'*H, eps1*E' + eps1*E1'; P3'*H, zeros(n,p)];
                H'*P2, H'*P3, -eps1 * eye(p), zeros(p,p);
                eps1*E + eps1*E1, zeros(p,n), zeros(p,p), -eps1*eye(p)] <= 0;
        lmi2 = [R, A1'*P2, A1'*P3, zeros(n,p), eps2*E1';
                P2'*A1, Z1, Z2, P2'*H, zeros(n,p);
                P3'*A1, Z2', Z3, P3'*H, zeros(n,p);
                zeros(p,n), H'*P2, H'*P3, eps2*eye(p), zeros(p);
                eps2*E1, zeros(p,n), zeros(p,n), zeros(p), eps2*eye(p)];
        
        const = [P1 >= small_scalar*eye(n), R >= small_scalar*eye(n), lmi1, lmi2];

        % Solve using the bisection method
        diagnostics = optimize(const,[],sdpsettings('solver','mosek','verbose',0));
        disp(diagnostics.problem)

        % Save if problem is feasible or not
        if diagnostics.problem == 0
            is_feasible_mat(i,j) = 1;
        end
    end
end

end

