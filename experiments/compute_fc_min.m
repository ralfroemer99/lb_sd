% Compute the minimum control frequency for an uncertain linear model
% learned from real-world measurements.

%% Get the data
load_path = "./experiments/data/from_lukas/";

A = readmatrix(load_path + "A.csv");
B = readmatrix(load_path + "B.csv");
% B = [-1;1];

Au = zeros(size(A));
Bu = zeros(size(B));

% Transform the matrix-valued uncertainty
[H,E,F] = transform_uncertainty(Au,Bu);

% Maximize the sampling time, i.e., minimize the control frequency.
eps_vec = logspace(-3,3,20);
% [Ts_max, K] = max_Ts_nominal(A,B);
[Ts_max, K, Ts_vec] = max_Ts_norm_bounded(A,B,H,E,F,eps_vec);

