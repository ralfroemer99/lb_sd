% Define the nonlinear model
switch system
    case 'quad1D'
        n = 3; m = 1;
        f = @quad1D_dynamics; params = quad1D_params;
        x_s = [0,0,0]; u_s = 0;
        sigma_n = [0.1,0.1,0.1];
        N_vec = 50:25:200;
        J1 = diag([100, 1, 100]); J2 = 0.1*eye(m);
    case 'quad2D'
        n = 6; m = 2;               % State and input dimensions
        f = @quad2D_dynamics; params = quad2D_params;
        x_s = [1,0,0,0,0,0]; u_s = [0.981,0.981];
        sigma_n = [0.1,0.1,0.1,0.1,0.1,0.1];
        N_vec = 200:200:600;       
        for_which_optimize = 1:9;
        J1 = diag([100,1,100,1,100,1]); J2 = 0.01*eye(m);

        delta_x0 = [0.2,0,0.2,0,0,0;
                    0.05,0,0.25,0,0,0;
                    0.25,0,0.05,0,0,0;
                    0.2,0,0.1,0,0.05,0;
                    0.05,0,0.2,0,-0.05,0];
    otherwise
        error('Wrong system type!')
end