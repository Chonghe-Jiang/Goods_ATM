%%% EC 2025 Submission
%%% ! Version: EC - Quasi Linear - Synthetic Data - Iterative method
%%% ! Optimal Version for Submission
%%% ! In QL Setting -> we do not do normalization anymore
clc
clear

% Define problem parameters
n = 50;  % Number of rows 
m = 50;   % Number of columns
B = ones(n, 1);  % Random B vector

% Todo: Change the threhold
max_iter = 1300;
max_iter_adaptive = 4500;
epsilon = 1e-4; 
plot_flag = true;

% Define the folder name
dataset_folder = 'synthetica_dataset_quasi';  % Changed folder name to include "quasi"
% Generate the filename for v based on n and m
%%% Todo: Remember to change the folder name if necessary
v_filename = sprintf('v_quasi_exp_%d_%d.mat', n, m);  % Changed filename to include "quasi"
v_filepath = fullfile(dataset_folder, v_filename); % Full path to the file
% Generate the filename for solver results based on n and m
solver_filename = sprintf('solver_quasi_exp_%d_%d.mat', n, m);  % Changed filename to include "quasi"
solver_filepath = fullfile(dataset_folder, solver_filename); % Full path to the file

% Check if the file exists. If it does, load 'v' from the file. Otherwise, generate 'v' and save it.
if exist(v_filepath, 'file') == 2
    load(v_filepath, 'v');  % Load 'v' from the file
    disp(['Loaded v from ', v_filepath]);
else
    % v = 20*rand(n, m)+5; 
    % v = randi([1, 10], n, m); % Draw valuations from integer uniform distribution
    v = exprnd(10, n, m); 
    % v = lognrnd(0, 10, n, m); % Draw valuations from log-normal distribution
    v = v ./ sum(v, 2); 
    save(v_filepath, 'v');  % Save 'v' to a file for future use
    disp(['Generated v and saved to ', v_filepath]);
end
% Check if the solver results file exists. If it does, load the results. Otherwise, solve and save the results.
if exist(solver_filepath, 'file') == 2
    load(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Load solver results from the file
    disp(['Loaded solver results from ', solver_filepath]);
else
    % Solve the problem using the solver
    [p_opt_solver, beta_opt, fval_solver, solve_time] = quasi_dual_solver(n, m, B, v);
    save(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Save solver results to a file
    disp(['Solved and saved solver results to ', solver_filepath]);
end
disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Transfer to a row vector

%%% ! Box constraint - Quasi Linear Version
p_lower = max(v .* B ./ (sum(abs(v),2)+B)); %%%!  quasi linear version
p_upper = max(v); %%%!  quasi linear version
mu_lower = log(p_lower);
mu_upper = log(p_upper);

%%% * Parameter for convexity and smoothness and stepsize
delta = 0.1;  
sigma = min(exp(mu_lower));
L = exp(max(mu_upper)) + (sum(B) / delta); 
adaptive = false;
step_size = 1e-3; %%% Todo: Subgradient stepsize
eta = 1;  %%% Todo: MD stepsize - now is the optimal step size 1

%%% * - ini of p0 and mu0
p0 = quasi_init_gd(p_lower,p_upper,sum(B));
% p0 = 1/2 * (p_upper + p_lower);
mu0 = log(p0); %
x0 = quasi_init_md(p0,B);

%%% * - solve the problem by subgradient
[solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = quasi_dual_subgradient(v, B, p0, max_iter, step_size, epsilon, plot_flag, p_opt_solver, fval_solver);
% disp(['Sub final gap: ', num2str(obj_values_sub(end))]);
disp(['Subgradient time: ', num2str(time_sub), ' seconds']);
disp(['Subgradient iterations: ', num2str(iter_sub)]);

%%% * - solve the problem by mirror descent
[solution_md, time_md, iter_md, obj_values_md, distance_md] = quasi_primal_md(v, B, x0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver);
disp(['MD iterations: ', num2str(iter_md)]);
disp(['MD time: ', num2str(time_md), ' seconds']);

%%% * - solve the problem by multiple-agd
adaptive_plot_flag = true;  % Set to true if you want to plot the results
plot_flag = false;
plot_flag_smooth = false;
adaptive = true;
phase_num = 20; %%% Todo: SGR - parameter
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive] = quasi_dual_adaptive(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive, phase_num);
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);
disp(['Adaptive AGD time: ', num2str(total_time_adaptive), ' seconds']);

%%% * - plot the descent graph
% Determine the smallest length among the three data sets
min_length = min([length(obj_values_md), length(obj_values_sub), length(obj_values_adaptive)]);

% Create x-axis values
x_md = 1:length(obj_values_md);
x_subgrad = 1:length(obj_values_sub);
x_adaptive = 1:length(obj_values_adaptive);

% Plotting
figure;
semilogy(x_md, abs(obj_values_md), '-d', 'DisplayName', 'Mirror Descent', 'LineWidth', 2); % Plot obj_mirror
hold on;
semilogy(x_subgrad, abs(obj_values_sub), '-d', 'DisplayName', 'Tatonnement', 'LineWidth', 2); % Plot obj_subgrad
semilogy(x_adaptive, abs(obj_values_adaptive), '-d', 'DisplayName', 'ATM', 'LineWidth', 2); % Plot obj_adaptive
hold off;

% Set font sizes and other properties
set(gca, 'FontSize', 15); % Set axis font size
xlabel('Iteration', 'FontSize', 25); % X-axis label with larger font size
ylabel('Objective Value Gap', 'FontSize', 25); % Y-axis label with larger font size
title(['n=', num2str(n), ', m=', num2str(m)], 'FontSize', 25); % Graph title with larger font size
legend show; % Show legend
grid on; % Enable grid

distance_md_to_solver = norm(solution_md - p_opt_solver); % Distance for Mirror Descent
distance_sub_to_solver = norm(solution_sub - p_opt_solver); % Distance for Subgradient
distance_adaptive_to_solver = norm(exp(solution_adaptive) - p_opt_solver); % Distance for Adaptive AGD

% Print the distances
disp(['Distance between Mirror Descent solution and solver solution: ', num2str(distance_md_to_solver)]);
disp(['Distance between Subgradient solution and solver solution: ', num2str(distance_sub_to_solver)]);
disp(['Distance between Adaptive AGD solution and solver solution: ', num2str(distance_adaptive_to_solver)]);