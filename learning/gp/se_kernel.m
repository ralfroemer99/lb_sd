function val = se_kernel(z1,z2,sigma_f,L_inv)
% Ensure that z1 and z2 are column vectors
if length(z1) ~= 1
    if size(z1,1) == 1
        z1 = z1';
    end
end
if length(z2) ~= 1
    if size(z2,1) == 1
        z2 = z2';
    end
end
% Compute the value of the squared exponential function
val = sigma_f^2 * exp(-1/2 * (z1 - z2)' * L_inv * (z1 - z2));
end

