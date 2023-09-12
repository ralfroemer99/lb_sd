function [Au,Bu] = dSigma2AB(dSigma,n)
% Transform the GP variance derivative into the uncertainty of a linearized 
% system such that it captures the true linearized dynamics with probability of at least (p_tilde)^n.

% Probability threshold per dimension
p_tilde = 0.99;

% Get input dimension
m = size(dSigma,2) - n;

C = zeros(n,n+m);
gam = sqrt(chi2inv(p_tilde,size(dSigma,2)));
for i = 1:n
    C(i,:) = gam*sqrt(diag(dSigma(:,:,i)));
end

Au = C(:,1:n);
Bu = C(:,n+1:end);
end

