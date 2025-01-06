% Test Program - Simple testing for the iterative method
clc
clear
% Set random seed for reproducibility
rng(42);

% Define problem parameters
n = 500;  % Number of rows
m = 500;   % Number of columns
B = rand(n, 1);  % Random B vector
v = rand(n, m);  % Random v matrix

% Set common parameters
max_iter = 5000;
epsilon = 1e-1; % Stopping criteria with epsilon
plot_flag = true;

% Set specific parameters
% AGD parameter
delta = 0.05;    
mu_lower = log(max(v .* B ./ sum(abs(v))));
mu_upper = log(norm(B, 1) * ones(1, m));
mu0 = (1/2)*(mu_lower+mu_upper); %
sigma = min(exp(mu_lower));
L = (exp(max(mu_upper)) + sum(B) / delta); 
adaptive = false;
% Subgradient parameters
p0 = exp(mu0);  % Initial point - Need to have the same for the later usage
step_size = 0.0001;
% Mirror descent algorithm
x_0 = rand(n, m);  % Initial point
eta = 0.1;  % Step size of the MD method

%%% Testing
% 1. Test linear_dual_solver
disp('Testing linear_dual_solver:');
[p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
disp(['Solver optimal value: ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Transfer to a row vector

% 2. Test linear_dual_subgradient
disp('Testing linear_dual_subgradient:');
[solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = linear_dual_subgradient(v, B, p0, max_iter, step_size, epsilon, plot_flag, p_opt_solver, fval_solver);
disp(['Sub final gap: ', num2str(obj_values_sub(end))]);
disp(['Sub time: ', num2str(time_sub), ' seconds']);
disp(['Sub iterations: ', num2str(iter_sub)]);
%%% Todo - stopping criteria of the subgradient method

% 3. Test linear_primal_md
% disp('Testing linear_primal_md:');
% [solution_md, time_md, iter_md, obj_values_md, distance_md] = linear_primal_md(v, B, x_0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver);
% disp(['MD final gap: ', num2str(obj_values_md(end))]);
% disp(['MD time: ', num2str(time_md), ' seconds']);
% disp(['MD iterations: ', num2str(iter_md)]);
%%% Todo - stepsize choice of the mirror descent method

% 4. Test linear_dual_agd
disp('Testing linear_dual_agd:');
[solution_agd, time_agd, iter_agd, obj_values_agd, dis_agd, convergence_agd] = linear_dual_agd(v, B, mu0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, p_opt_solver, fval_solver, adaptive);
disp(['AGD final gap: ', num2str(obj_values_agd(end))]);
disp(['AGD time: ', num2str(time_agd), ' seconds']);
disp(['AGD iterations: ', num2str(iter_agd)]);
%%% Todo - parameter choice in the accelerated method

% 5. Test linear_dual_adaptive
disp('Testing linear_dual_adaptive:');
adaptive_plot_flag = true;  % Set to true if you want to plot the results
plot_flag = false;
adaptive = true;
phase_num = 5;
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_adaptive, dis_adaptive] = linear_dual_adaptive(v, B, mu0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, p_opt_solver, fval_solver, adaptive, phase_num);
disp(['Adaptive AGD final gap: ', num2str(obj_adaptive(end))]);
disp(['Adaptive AGD total time: ', num2str(total_time_adaptive), ' seconds']);
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);


%%% Compare results
%%% Todo - Compare what - time, function value gap, whatever can be the metric
% disp('Results comparison:');
% disp(['Solver optimal value: ', num2str(fval_solver)]);
% disp(['Sub final value gap: ', num2str(obj_values_sub(end))]);
% disp(['AGD final value gap: ', num2str(obj_values_agd(end))]);
% disp(['Adaptive AGD final value gap: ', num2str(obj_adaptive(end))]);
% disp(['Mirror descent final value gap: ', num2str(obj_values_md(end))]);