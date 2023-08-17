function Delta = Sigma2confidenceInterval(Sigma)
    % Compute a matrix-valued confidence interval based on the covariance
    % matrices of the normally distributed rows.
    % TODO: Consider the off-diagonal elements of the covariance matrix!
    % dim(Sigma{i}) = (n,n) or (n,m)
    gam = 3;
    [n,m] = size(Sigma{1});
    Delta = zeros(n,m);
    for i = 1:n
        Delta(i,:) = gam * sqrt(diag(Sigma{i}));
    end
end

