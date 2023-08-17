function [A,B] = dmu2AB(dmu,n)
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

