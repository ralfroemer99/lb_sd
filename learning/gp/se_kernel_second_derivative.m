function d2kdz2 = se_kernel_second_derivative(sigma_f,L_inv)
% Compute the second derivative of the squared-exponential kernel. 
% See Berkenkamp and Schoellig, 2015, eq. (9). 
    d2kdz2 = sigma_f^2 * L_inv;
end

