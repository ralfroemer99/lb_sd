function [sigma_f,L] = get_hyperparameters(gpr_model)
% Extract the hyperparameters of a gpr object.

d = size(gpr_model.X,2);
switch gpr_model.KernelFunction
    case 'SquaredExponential'
        params = gpr_model.KernelInformation.KernelParameters;
        L = diag(ones(d,1)*params(1));          % length scale
        sigma_f = params(2);
    case 'ARDSquaredExponential'
        params = gpr_model.KernelInformation.KernelParameters;
        L = diag(params(1:d));
        sigma_f = params(end);
    otherwise
        error('Kernel must be of squared exponential type!');
end
end

