function [x2,inBetween] = shaded_plot_min_max(x,y_min,y_max)
    % Helper function to draw plots with shaded +/- one standard deviation
    % Dimensions: dim(x) = [1,n]
    x2 = [x, fliplr(x)];
    inBetween = [y_max, fliplr(y_min)];
end

