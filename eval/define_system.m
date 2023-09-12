% Define the nonlinear model as a 1-D or 2-D quadrotor
switch system
    case 'quad1D'
        % State and input dimensions
        n = 3; m = 1;

        % System dynamics and parameters
        f = @quad1D_dynamics; params = quad1D_params;

        % Equilibrium/operating point
        x_s = [0,0,0]; u_s = 0;

        % Observation noise standard deviation
        sigma_n = [0.1,0.1,0.1];

        % Number of datapoints to train from
        N_vec = 50:25:200;

        % Cost weight matrices
        J1 = diag([100, 1, 100]); J2 = 0.1*eye(m);
    case 'quad2D'
        % State and input dimensions
        n = 6; m = 2;

        % System dynamics and parameters
        f = @quad2D_dynamics; params = quad2D_params;

        % Equilibrium/operating point
        x_s = [1,0,0,0,0,0]; u_s = [0.981,0.981];

        % Observation noise standard deviation
        sigma_n = [0.1,0.1,0.1,0.1,0.1,0.1];

        % Number of datapoints to train from
        N_vec = 200:200:800;       
        
        % Cost weight matrices
        J1 = diag([100,1,100,1,100,1]); J2 = 0.01*eye(m);
        
        % Initial conditions to simulate
        delta_x0 = [0.2,0,0.2,0,0,0;
                    0.05,0,0.25,0,0,0;
                    0.25,0,0.05,0,0,0;
                    0.2,0,0.1,0,0.05,0;
                    0.05,0,0.2,0,-0.05,0];
    otherwise
        error('Undefined system type!')
end