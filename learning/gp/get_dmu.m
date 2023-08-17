function dmu = get_dmu(sigma_f,L_inv,K,y,Z,z_star)
% Compute the derivative of the GP mean function at a test point
% TODO: Documentation
dmu = se_kernel_derivative(sigma_f,L_inv,Z,z_star) * (K \ y);
end

