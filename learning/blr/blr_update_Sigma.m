function Sigma_new = blr_update_Sigma(z,Sigma_old,sigma_n)
    % Update the covariance of the parameter estimate in BLR
    Sigma_new = inv((inv(Sigma_old) + 1/sigma_n^2 * (z * z')));
end

