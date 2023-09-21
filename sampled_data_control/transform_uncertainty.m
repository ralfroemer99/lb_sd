function [H, E, F] = transform_uncertainty(A_u, B_u)
% Get dimensions
[n, m] = size(B_u);

% Initialize
E = zeros(n^2 + n*m, n);
F = zeros(n^2 + n*m, m);

% Build H
H = [kron(eye(n), ones(1,n)), kron(eye(n), ones(1,m))];

% Build E
for i = 1:n
    E(1+(i-1)*n:i*n,:) = diag(A_u(i,:));
end

% Build E1
for i = 1:n
    F(n^2+1+(i-1)*m:n^2+i*m,:) = diag(B_u(i,:));
end

end

