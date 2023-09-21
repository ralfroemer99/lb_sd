function [idmodel] = identify_ss_model(data, delay, discrete, use_acc)
    %% Set up data
    timestamps = data.timestamps;
    states = data.states;
    inputs = data.inputs;
    dt = data.dt;
    
    state_dim = size(states, 1);
    input_dim = size(inputs, 1);

    % Set the initial state to zero    
    % states(1, :) = states(1, :) - states(1, 1);
    % states(2, :) = states(2, :) - states(2, 1);

    %% Provide discrete-time LTI State space system structure
    output_dim = state_dim;
    
    % Set dummy matrices    
    B = zeros(state_dim, input_dim); 
    C = eye(state_dim);
    D = zeros(output_dim, input_dim);
    
    % Initialize dummy state space system
    if discrete
        A = eye(state_dim, state_dim);
        idmodel = idss(A, B, C, D, 'Ts', dt);
        idmodel.Structure.A.Value(1, 1) = 1;
        idmodel.Structure.A.Free(1,1) = false;
        idmodel.Structure.A.Value(1, 2) = dt;
        idmodel.Structure.A.Free(1,2) = false;
        
        idmodel.Structure.A.Value(2, 1) = 0;
        idmodel.Structure.A.Free(2, 1) = false;
        idmodel.Structure.A.Value(2, 2) = 1;
        idmodel.Structure.A.Free(2, 2) = false;
        
        if use_acc
            idmodel.Structure.A.Value(1, 3) = dt * dt / 2;
            idmodel.Structure.A.Free(1, 3) = false;
            idmodel.Structure.A.Value(2, 3) = dt;
            idmodel.Structure.A.Free(2, 3) = false;
            idmodel.Structure.A.Value(3, 1) = 0;
            idmodel.Structure.A.Free(3, 1) = false;
            idmodel.Structure.A.Value(3, 2) = 0;
            idmodel.Structure.A.Free(3, 2) = false;
            idmodel.Structure.A.Value(3, 3) = dt;
            idmodel.Structure.A.Free(3, 3) = false;
            idmodel.Structure.A.Value(3, 3) = 0;
            idmodel.Structure.A.Free(3, 3) = false;
    
        %     idmodel.Structure.B.Value(3, 1) = 18;
        %     idmodel.Structure.B.Free(3, 1) = false;
        % else
        %     idmodel.Structure.B.Value(2, 1) = 18;
        %     idmodel.Structure.B.Free(2, 1) = false;
        end
    else
        A = zeros(state_dim, state_dim);
        idmodel = idss(A, B, C, D, 'Ts', 0);
        idmodel.Structure.A.Value(1, 1) = 0;
        idmodel.Structure.A.Free(1,1) = false;
        idmodel.Structure.A.Value(1, 2) = 1;
        idmodel.Structure.A.Free(1,2) = false;
        idmodel.Structure.A.Value(2, 1) = -0.0933;
        idmodel.Structure.A.Free(2, 1) = false;
        idmodel.Structure.A.Value(2, 2) = 0.1035;
        idmodel.Structure.A.Free(2, 2) = false;

        if use_acc
            idmodel.Structure.A.Value(1, 3) = 0;
            idmodel.Structure.A.Free(1, 3) = false;
            idmodel.Structure.A.Value(2, 3) = 1;
            idmodel.Structure.A.Free(2, 3) = false;
            idmodel.Structure.A.Value(3, 1) = 0;
            idmodel.Structure.A.Free(3, 1) = false;
            idmodel.Structure.A.Value(3, 2) = 0;
            idmodel.Structure.A.Free(3, 2) = false;
            idmodel.Structure.A.Value(3, 3) = 0;
            idmodel.Structure.A.Free(3, 3) = false;
            idmodel.Structure.A.Value(3, 3) =  1;
            idmodel.Structure.A.Free(3, 3) = false;
    
            % idmodel.Structure.B.Value(3, 1) = 18 / dt;
            % idmodel.Structure.B.Free(3, 1) = false;
        else
            idmodel.Structure.B.Value(2, 1) = 18.0919;
            idmodel.Structure.B.Free(2, 1) = false;
        end
    end

    % Set constraints on the system matrices:
    % - full state measurement from the simulator --> C is identity
    % - D is zero
    % - the first entry of B is zero, because inputs do not directly affect the
    %   position
    for i = 1 : state_dim
        for j = 1 : state_dim
            if i == j
                idmodel.Structure.C.Value(i, j) = 1;
            else
                issestdmodel.Structure.C.Value(i, j) = 0;
            end
            idmodel.Structure.C.Free(i, j) = false;
        end
        for j = 1 : input_dim
            idmodel.Structure.D.Free(i,j) = false;
        end
    end

    idmodel.Structure.B.Value(1, 1) = 0;
    idmodel.Structure.B.Free(1,1) = false;
    if use_acc
        idmodel.Structure.B.Value(2, 1) = 0;
        idmodel.Structure.B.Free(2, 1) = false;
    end

    %% Identify discrete-time LTI system
    % id_data = iddata(states', inputs', dt);
    % idmodel = ssest(id_data, idmodel, 'InputDelay', delay);
    id_data = iddata(states', inputs', dt);
    idmodel = ssest(id_data, idmodel);
    A = idmodel.A;
    B = idmodel.B;
    C = idmodel.C;

    %% Simulate model using control input data
    verify_model(idmodel, data, discrete);
end