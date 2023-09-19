# Learning-based sampled-data control
This repository provides the code to reproduce the results in the paper "Is Data All That Matters? The Role of Control Frequency for the Stability and Closed-Loop Performance of Uncertain Systems", available HERE.

In this paper, we combine Gaussian process (GP) model learning with sampled-data control to study the role of the control frequency for uncertain systems. Our proposed framework has the following structure:
![alt text](https://github.com/ralfroemer99/lb_sd/blob/main/framework_structure.jpg?raw=true)

## Installation
- Install MATLAB with the Statistics and Machine Learning Toolbox (https://matlab.mathworks.com/).
- Clone the repository.
- Install the MOSEK Optimization suite (https://docs.mosek.com/latest/install/installation.html).
- Download YALMIP (https://yalmip.github.io/download/).

## Run the Code
### Compute the Minimum Control Frequency (MCF)
- Run eval/min_fc_over_N.m


The script creates datasets with different amounts of randomly drawn samples. For each dataset, the dynamics are learned using GP regression and linearized. For each linearized system, the MCF is computed by solving a semidefinite program.

The resulting plot corresponds to Figure 2 of the paper.
![alt text](https://github.com/ralfroemer99/lb_sd/blob/main/fig_2.jpg?raw=true)


### Simulate the System at Different Control Frequencies
- Step 1: Run eval/sim_trajectories_step1.m
- Step 2: Run eval/sim_trajectories_step2.m


The first script creates datasets with different amounts of randomly drawn samples and learns and linearizes the dynamics similar to eval/min_fc_over_N.m. For each linearized system, the MCF is computed by solving a semidefinite program. Then, optimized controllers operating at {1.25,1.5,2}*MCF are computed by solving another semidefinite program. The second script simulates the 2D quadrotor with each controller and the corresponding control frequency for different amounts of training data.

The resulting plot corresponds to Figure 3 of the paper.
![alt text](https://github.com/ralfroemer99/lb_sd/blob/main/fig_3.jpg?raw=true)

### Analyze the Tradeoff between Data and Control Frequency
- Step 1: Run eval/analyze_performance_step1.m
- Step 2: Run eval/analyze_performance_step2.m


The first script creates datasets with different amounts of randomly drawn samples and learns and linearizes the dynamics similar to eval/min_fc_over_N.m. Optimized controllers operating at fixed frequencies between 10 and 30 Hz are computed by solving a semidefinite program. The second script simulates the 2D quadrotor system with each controller and control frequency for different initial conditions and computes the average cost.

The resulting plot corresponds to Figure 4 of the paper.
![alt text](https://github.com/ralfroemer99/lb_sd/blob/main/fig_4.jpg?raw=true)
