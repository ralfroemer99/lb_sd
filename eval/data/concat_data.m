% Concatenate the data obtained from different simulations and seeds
clear

% Provide load and save path
load_path_base = "./eval/data/analyze_performance/main_fig_low_res/";
save_path = "./eval/data/analyze_performance/main_fig_low_res/all/";

% Provide list of subpaths
% paths = {'1to5_seed0','6_seed1','7to8_seed2','9_seed3'};
paths = {"1to5_seed0","1to5_seed0","1to5_seed0"};

% Initialize variables
load(load_path_base + paths{1} + "/" + "A_all",'A_all');
load(load_path_base + paths{1} + "/" + "B_all",'B_all');
load(load_path_base + paths{1} + "/" + "Au_all",'Au_all');
load(load_path_base + paths{1} + "/" + "Bu_all",'Bu_all');
load(load_path_base + paths{1} + "/" + "N_vec",'N_vec');
load(load_path_base + paths{1} + "/" + "K_all",'K_all');
load(load_path_base + paths{1} + "/" + "fc_vec",'fc_vec');

for i = 2:length(paths)
    % Load next variables
    A_all_tmp = load(load_path_base + paths{i} + "/" + "A_all").A_all;
    B_all_tmp = load(load_path_base + paths{i} + "/" + "B_all",'B_all').B_all;
    Au_all_tmp = load(load_path_base + paths{i} + "/" + "Au_all",'Au_all').Au_all;
    Bu_all_tmp = load(load_path_base + paths{i} + "/" + "Bu_all",'Bu_all').Bu_all;
    K_all_tmp = load(load_path_base + paths{i} + "/" + "K_all",'K_all').K_all;
    
    % Concatenate
    A_all = cat(4,A_all,A_all_tmp);
    B_all = cat(4,B_all,B_all_tmp);
    Au_all = cat(4,Au_all,Au_all_tmp);
    Bu_all = cat(4,Bu_all,Bu_all_tmp);
    K_all = cat(4,K_all,K_all_tmp);
end

% Save variables
save(save_path + "A_all",'A_all');
save(save_path + "B_all",'B_all');
save(save_path + "Au_all",'Au_all');
save(save_path + "Bu_all",'Bu_all');
save(save_path + "N_vec",'N_vec');
save(save_path + "K_all",'K_all');
save(save_path + "fc_vec",'fc_vec');

