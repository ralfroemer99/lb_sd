function K = get_K(Z,sigma_f,L_inv,sigma_n)
% Compute the kernel gram matrix for the inputs in Z and the 
% squared-exponential kernel.

[N,~] = size(Z);
    K = zeros(N);
    for i = 1:N
        z_i = Z(i,:);
        for j=1:N
            z_j = Z(j,:);
            K(i,j) = se_kernel(z_i,z_j,sigma_f,L_inv);
        end
    end
    K = K + sigma_n^2 * eye(N);
end