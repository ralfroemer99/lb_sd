clear
close all

% Define real system
A = [1 1; 2 -1];
B = [0; 1];

% Dimensions
[n,m] = size(B);

% Define control frequency
Ts1 = 0.05:0.05:1;
alpha1 = zeros(1,length(Ts1));

% Compute maximum bound on the uncertainty
for i = 1:length(Ts1)
    alpha1(i) = max_disturbance(A,B,sqrt(n)*eye(n),Ts1(i));
end

% alpha2 = 0.05:0.05:1;
% Ts2 = zeros(1,length(alpha2));
% for i = 1:length(alpha2)
%     Ts2(i) = max_sampling_time(A,B,sqrt(n)*eye(n),alpha2(i));
% end

plot(Ts1,alpha1,'r'); hold on
% plot(Ts2,alpha2,'g')
