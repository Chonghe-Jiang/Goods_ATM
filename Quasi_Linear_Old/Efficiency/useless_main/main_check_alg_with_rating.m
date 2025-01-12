% Test Program - Simple testing for the iterative method
%%% Todo - The key lies in the same initialization
clc
clear

% Set the working directory to the folder containing the Dataset folder
cd('C:\Users\s1155203585\Desktop\EG_EXP\Linear\v0_funding_supp\Dataset');
% Read the CSV file into a MATLAB matrix
csvFilePath = fullfile('Dataset', 'Ratings.csv');
v = readmatrix(csvFilePath);
v = v/norm(v,"fro");
[n,m] = size(v);
B = ones(n,1);
B = B/norm(B,"fro");

% Now you can use dataMatrix in your MATLAB code

% Set common parameters
max_iter = 10000;
max_iter_adaptive = 4500;
epsilon = 1e-1; % Stopping criteria with epsilon
plot_flag = true;

%%% Todo - Box constraint  
p_lower = max(v .* B ./ sum(abs(v)));
p_upper = norm(B, 1) * ones(1, m);
mu_lower = log(max(v .* B ./ sum(abs(v))));
mu_upper = log(norm(B, 1) * ones(1, m));

%%% Todo - Parameter for convexity and smoothness and stepsize
delta = 0.3;  
sigma = min(exp(mu_lower));
L = (exp(max(mu_upper)) + sum(B) / delta); 
adaptive = false;
step_size = 0.0001;
eta = 0.1;  

%%% Todo - ini of p0 and mu0
p0 = linear_init_gd(p_lower,p_upper,sum(B));
mu0 = log(p0); %
x0 = linear_init_md(p0,B);

%%% Testing
% 1. Test linear_dual_solver
% disp('Testing linear_dual_solver:');
[p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
% disp(['Solver optimal value: ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Transfer to a row vector

% 2. Test linear_dual_subgradient
% disp('Testing linear_dual_subgradient:');
% [solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = linear_dual_subgradient(v, B, p0, max_iter, step_size, epsilon, plot_flag, p_opt_solver, fval_solver);
% disp(['Sub final gap: ', num2str(obj_values_sub(end))]);
% disp(['Sub time: ', num2str(time_sub), ' seconds']);
% disp(['Sub iterations: ', num2str(iter_sub)]);
%%% Todo - stopping criteria of the subgradient method

% 3. Test linear_primal_md
% disp('Testing linear_primal_md:');
[solution_md, time_md, iter_md, obj_values_md, distance_md] = linear_primal_md(v, B, x0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver);
% disp(['MD final gap: ', num2str(obj_values_md(end))]);
disp(['MD iterations: ', num2str(iter_md)]);
disp(['MD time: ', num2str(time_md), ' seconds']);

% 4. Test linear_dual_agd
% disp('Testing linear_dual_agd:');
% [solution_agd, time_agd, iter_agd, obj_values_agd, dis_agd, convergence_agd] = linear_dual_agd(v, B, mu0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, p_opt_solver, fval_solver, adaptive);
% % disp(['AGD final gap: ', num2str(obj_values_agd(end))]);
% disp(['AGD iterations: ', num2str(iter_agd)]);
% disp(['AGD time: ', num2str(time_agd), ' seconds']);

%%% Todo - parameter choice in the accelerated method

% 5. Test linear_dual_adaptive
% disp('Testing linear_dual_adaptive:');
adaptive_plot_flag = true;  % Set to true if you want to plot the results
plot_flag = false;
adaptive = true;
phase_num = 5;
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_adaptive, dis_adaptive] = linear_dual_adaptive(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, p_opt_solver, fval_solver, adaptive, phase_num);
% disp(['Adaptive AGD final gap: ', num2str(obj_adaptive(end))]);
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);
disp(['Adaptive AGD time: ', num2str(total_time_adaptive), ' seconds']);



%%% Compare results
%%% Todo - Compare what - time, function value gap, whatever can be the metric
% disp('Results comparison:');
% disp(['Solver optimal value: ', num2str(fval_solver)]);
% disp(['Sub final value gap: ', num2str(obj_values_sub(end))]);
% disp(['AGD final value gap: ', num2str(obj_values_agd(end))]);
% disp(['Adaptive AGD final value gap: ', num2str(obj_adaptive(end))]);
% disp(['Mirror descent final value gap: ', num2str(obj_values_md(end))]);