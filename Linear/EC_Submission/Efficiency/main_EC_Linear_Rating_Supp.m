%%% EC 2025 Submission
%%% ! Version: EC - Linear - Rating dataset - Iterative method - Unit budget (every element in B is 1)
%%% ! Optimal Version for Submission
clc
clear

% Load dataset
v = readmatrix('Dataset/Ratings_kroer.csv') + 0.1;
v = floor(v) + 1; % Preprocess the data
[n, m] = size(v);
B = ones(n, 1); % Unit budget

% Set common parameters
max_iter = 20000;
max_iter_adaptive = 4500;
epsilon = 1e-4; % Stopping criteria with epsilon
plot_flag = true;

%%% * Box constraint  
p_lower = max(v .* B ./ sum(abs(v), 2)); 
p_upper = norm(B, 1) * ones(1, m);
mu_lower = log(p_lower);
mu_upper = log(p_upper);

%%% * Parameter for convexity and smoothness and stepsize
delta = 0.1;  % Smoothing parameter
sigma = min(exp(mu_lower));
L = exp(max(mu_upper)) + (sum(B) / delta); 
adaptive = false;
step_size = 1e-4; % Subgradient stepsize
eta = 0.5;  % MD stepsize

%%% * - Initialize p0 and mu0
p0 = linear_init_gd(p_lower, p_upper, sum(B));
mu0 = log(p0); 
x0 = linear_init_md(p0, B);

%%% * Load solver results
load('/Users/chjiang/Dropbox/EG_EXP/Linear/EC_Submission/Efficiency/Solver_Rating.mat', 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');

%%% * Solve the problem using Subgradient
[solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = linear_dual_subgradient(v, B, p0, max_iter, step_size, epsilon, plot_flag, p_opt_solver, fval_solver);
disp(['Subgradient time: ', num2str(time_sub), ' seconds']);
disp(['Subgradient iterations: ', num2str(iter_sub)]);

%%% * Solve the problem using Adaptive Subgradient
switch_step_sub = 5000; % Number of iterations after which the step size decreases
[solution_sub_adaptive, obj_values_sub_adaptive, dis_sub_adaptive, time_sub_adaptive, iter_sub_adaptive] = linear_dual_subgradient_adaptive(v, B, p0, max_iter, step_size, epsilon, plot_flag, p_opt_solver, fval_solver, switch_step_sub);
disp(['Adaptive Subgradient time: ', num2str(time_sub_adaptive), ' seconds']);
disp(['Adaptive Subgradient iterations: ', num2str(iter_sub_adaptive)]);

%%% * Solve the problem using Mirror Descent
[solution_md, time_md, iter_md, obj_values_md, distance_md] = linear_primal_md(v, B, x0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver);
disp(['MD iterations: ', num2str(iter_md)]);
disp(['MD time: ', num2str(time_md), ' seconds']);

%%% * Solve the problem using Adaptive Mirror Descent
switch_step_md = 5000; % Number of iterations after which the step size decreases
[solution_md_adaptive, time_md_adaptive, iter_md_adaptive, obj_values_md_adaptive, distance_md_adaptive] = linear_primal_md_adaptive(v, B, x0, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver, switch_step_md);
disp(['Adaptive MD iterations: ', num2str(iter_md_adaptive)]);
disp(['Adaptive MD time: ', num2str(time_md_adaptive), ' seconds']);

%%% * Solve the problem using Adaptive AGD
adaptive_plot_flag = true;  
plot_flag = false;
plot_flag_smooth = false;
adaptive = true;
phase_num = 30;
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive] = linear_dual_adaptive(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive, phase_num);
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);
disp(['Adaptive AGD time: ', num2str(total_time_adaptive), ' seconds']);

%%% * Plot the descent graph
x_md = 1:length(obj_values_md);
x_subgrad = 1:length(obj_values_sub);
x_subgrad_adaptive = 1:length(obj_values_sub_adaptive);
x_adaptive = 1:length(obj_values_adaptive);
x_md_adaptive = 1:length(obj_values_md_adaptive);

figure;
semilogy(x_subgrad, abs(obj_values_sub), '-d', 'DisplayName', 'Tatonnement', 'LineWidth', 2); % Plot Subgradient
hold on;
semilogy(x_subgrad_adaptive, abs(obj_values_sub_adaptive), '-d', 'DisplayName', 'Adaptive Tatonnement', 'LineWidth', 2); % Plot Adaptive Subgradient
semilogy(x_md, abs(obj_values_md), '-d', 'DisplayName', 'Mirror Descent', 'LineWidth', 2); % Plot Mirror Descent
semilogy(x_md_adaptive, abs(obj_values_md_adaptive), '-d', 'DisplayName', 'Adaptive MD', 'LineWidth', 2); % Plot Adaptive MD
semilogy(x_adaptive, abs(obj_values_adaptive), '-d', 'DisplayName', 'SGA', 'LineWidth', 2); % Plot Adaptive AGD
hold off;

% Set font sizes and other properties
set(gca, 'FontSize', 15); % Set axis font size
xlabel('Iteration', 'FontSize', 25); % X-axis label with larger font size
ylabel('Objective Value Gap', 'FontSize', 25); % Y-axis label with larger font size
title('Movie Rating Dataset', 'FontSize', 25); % Graph title with larger font size
legend('show', 'Location', 'best'); % Show legend and position it at the best location
grid on; % Enable grid

%%% * Calculate distances to solver solution
distance_md_to_solver = norm(solution_md - p_opt_solver); % Distance for Mirror Descent
distance_sub_to_solver = norm(solution_sub - p_opt_solver); % Distance for Subgradient
distance_sub_adaptive_to_solver = norm(solution_sub_adaptive - p_opt_solver); % Distance for Adaptive Subgradient
distance_adaptive_to_solver = norm(exp(solution_adaptive) - p_opt_solver); % Distance for Adaptive AGD
distance_md_adaptive_to_solver = norm(solution_md_adaptive - p_opt_solver); % Distance for Adaptive MD

% Print the distances
disp(['Distance between Mirror Descent solution and solver solution: ', num2str(distance_md_to_solver)]);
disp(['Distance between Subgradient solution and solver solution: ', num2str(distance_sub_to_solver)]);
disp(['Distance between Adaptive Subgradient solution and solver solution: ', num2str(distance_sub_adaptive_to_solver)]);
disp(['Distance between Adaptive AGD solution and solver solution: ', num2str(distance_adaptive_to_solver)]);
disp(['Distance between Adaptive MD solution and solver solution: ', num2str(distance_md_adaptive_to_solver)]);