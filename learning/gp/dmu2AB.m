function [A,B] = dmu2AB(dmu,n)
% Transform the GP mean derivative into the system matrices of the nominal
% linearized system.

if size(dmu,2) ~= n
    dmu = dmu';
end
if size(dmu,2) ~= n
    error('Wrong dimensions!');
end
m = size(dmu,1) - n;
tmp = zeros(n,n+m);
for i = 1:n
    tmp(i,:) = dmu(:,i)';
end
A = tmp(:,1:n);
B = tmp(:,n+1:n+m);
end

