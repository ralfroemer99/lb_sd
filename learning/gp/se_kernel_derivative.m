function dkdz = se_kernel_derivative(sigma_f,L_inv,Z,z_star)
% Compute the derivative of the squared-exponential kernel. 
% See Berkenkamp and Schoellig, 2015, eq. (8). 

% Ensure that z_star is a column vector
z_star = z_star(:);

[N,n] = size(Z);
mat = zeros(n,N);
for i = 1:N
    mat(:,i) = (Z(i,:)'-z_star) * se_kernel(Z(i,:),z_star,sigma_f,L_inv);
end
dkdz = L_inv*mat;
end

