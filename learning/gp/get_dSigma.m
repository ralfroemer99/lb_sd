function dSigma = get_dSigma(sigma_f,L_inv,K,Z,z_star)
% Compute the derivative of the GP variance at a test point z_star. 

dkdz = se_kernel_derivative(sigma_f,L_inv,Z,z_star);
dSigma = se_kernel_second_derivative(sigma_f,L_inv) - dkdz * (K \ dkdz');
end

