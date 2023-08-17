function mu_new = blr_update_mu(z,y,mu_old,Sigma_old,Sigma_new,sigma_n)
    % Update mean of the posterior parameter distribution of BLR
    mu_new = Sigma_new * (Sigma_old\mu_old + 1/sigma_n^2 * z * y);
end

