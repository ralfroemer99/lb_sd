function [A_new, Sigma_new] = blr_update_A(x,u,dx,A_old,B,Sigma_old,sigma_n)
    % Get dimensions and initialize variables
    n = size(A_old,1);
    A_new = zeros(n);
    Sigma_new = cell(1,n);

    % Get training input and target
    z = x;
    y = dx - B * u;

    % Update the mean and the covariance of the rows of A
    for i = 1:n
        % Update variance
        Sigma_new{i} = blr_update_Sigma(x,Sigma_old{i},sigma_n);
        A_new(i,:) = blr_update_mu(z,y(i),A_old(i,:)',Sigma_old{i},Sigma_new{i},sigma_n)';
    end
end

