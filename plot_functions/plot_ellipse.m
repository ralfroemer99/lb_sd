function plot_ellipse(P, mu)
 % plots an ellipse of the form xEx = 1
 R = chol(P);
 t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
 z = [cos(t); sin(t)];
 ellipse = inv(R) * z;
 ellipse(1,:) = ellipse(1,:) + mu(1);
 ellipse(2,:) = ellipse(2,:) + mu(2);
 plot(ellipse(1,:), ellipse(2,:))
end 

