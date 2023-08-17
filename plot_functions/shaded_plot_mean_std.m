function [x2,inBetween] = shaded_plot_mean_std(x,y_mean,y_std)
    % Helper function to draw plots with shaded +/- one standard deviation
    curve1 = y_mean + y_std;
    curve2 = y_mean - y_std;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
end

