# lb_sd
This repository provides the code for the results reported in the paper "The role of Control Frequency for Uncertain Systems" available HERE.

## Installation
- Install MATLAB with the Statistics and Machine Learning Toolbox
- Clone the repository
- Install the MOSEK Optimization suite (available at https://docs.mosek.com/latest/install/installation.html)
- Download YALMIP (available at https://yalmip.github.io/download/) 

## Run the Code
# Compute the minimum control frequency (MCF)
- Run eval/min_fc_over_N.m
The script creates datasets with different amounts of randomly drawn samples. For each dataset, the dynamics are learned using GP regression and linearized. For each linearized system, the MCF is computed by solving a semidefinite program.
The resulting plot corresponds to Fig. 2 of the paper.


# Simulate the system at different control frequencies
- (Optional:) Run eval/sim_trajectories_step1.m
- Run eval/sim_trajectories_step2.m
The first script creates datasets with different amounts of randomly drawn samples similar to eval/min_fc_over_N.m. For each dataset, the dynamics are learned using GP regression and linearized. For each linearized system, the MCF is computed by solving a semidefinite program. Then, optimized controllers operating at {1.25,1.5,2}*MCF are computed by solving another semidefinite program. The quadrotor system is simulated with each controller and the corresponding control frequency for different amounts of data.
The resulting plot corresponds to Fig. 3 of the paper.

# Analyze the tradeoff between data and control frequency in terms of performance
- (Optional:) Run eval/analyze_performance_step1.m
- Run eval/analyze_performance_step2.m
The first script creates datasets with different amounts of randomly drawn samples similar to eval/min_fc_over_N.m. For each dataset, the dynamics are learned using GP regression and linearized. Optimized controllers operating at fixed frequencies between 10 and 30 Hz are computed by solving a semidefinite program. The quadrotor system is simulated with the different controllers and control frequencies for different initial conditions, and the average cost is computed.
The resulting plot corresponds to Fig. 4 of the paper.
