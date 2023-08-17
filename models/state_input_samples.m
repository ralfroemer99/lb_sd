function [x_rand,u_rand] = state_input_samples(xlim_low,xlim_high,ulim_low,ulim_high,N)
    % dim(xlim_low) = dim(xlim_high) = (n,1)
    % dim(ulim_low) = dim(ulim_high) = (m,1)

    % Get dimensions
    n = length(xlim_low);
    m = length(ulim_low);

    % Create samples
    x_rand = zeros(n,N);
    u_rand = zeros(m,N);
    
    for i = 1:N
        % State
        for k = 1:n
            x_rand(k,i) = xlim_low(k) + rand * (xlim_high(k) - xlim_low(k));
        end
        % Input
        for k = 1:m
            u_rand(k,i) = ulim_low(k) + rand * (ulim_high(k) - ulim_low(k));
        end
    end
end

