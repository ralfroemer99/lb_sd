function [Au,Bu] = dSigma2AB(dSigma,n)
% dim(dSigma) = (n,n+m,n+m)
m = size(dSigma,2) - n;
tmp = zeros(n,n+m);

gam = sqrt(chi2inv(0.99,size(dSigma,2)));
for i = 1:n
    tmp(i,:) = gam*sqrt(diag(dSigma(:,:,i)));
end

Au = tmp(:,1:n);
Bu = tmp(:,n+1:end);
end

