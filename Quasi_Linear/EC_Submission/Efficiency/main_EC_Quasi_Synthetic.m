%%% EC 2025 Submission
%%% ! Version: EC - Linear - Sythetic Data - Iterative method
%%% ! Optimal Version for Submission
clc
clear
% Set random seed for reproducibility

%%% Todo: Start from the basic size
% Define problem parameters
n = 5;  % Number of rows 
m = 5;   % Number of columns
B = ones(n, 1);  % Random B vector
v = rand(n,m);
v = v ./ sum(v, 2); 

% Set common parameters
max_iter = 10000;
max_iter_adaptive = 4500;
epsilon = 1e-1; % Stopping criteria with epsilon %%% Todo: Change the threhold
plot_flag = true;

%%% ! Box constraint - Quasi Linear Version
p_lower = max(v .* B ./ (sum(abs(v),2)+B)); %%%!  quasi linear version
p_upper = max(max(v))*ones(1,m); %%%!  quasi linear version
mu_lower = log(p_lower);
mu_upper = log(p_upper);

%%% * Parameter for convexity and smoothness and stepsize
delta = 0.3;  
sigma = min(exp(mu_lower));
L = exp(max(mu_upper)) + (sum(B) / delta); 
adaptive = false;
step_size = 1e-5; %%% Todo: Subgradient stepsize
eta = 0.2;  %%% Todo: MD stepsize

%%% * - ini of p0 and mu0
p0 = quasi_init_gd(p_lower,p_upper,sum(B));
mu0 = log(p0); %
x0 = quasi_init_md(p0,B);

%%% * - solve the problem by solver
%%% ! No need for change
[p_opt_solver, beta_opt, fval_solver, solve_time] = quasi_dual_solver(n, m, B, v);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Transfer to a row vector

%%% ! Extra check: initialization
%%% * - solve the problem by subgradient
%%% ! Require changes for subgradient and function value
[solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = quasi_dual_subgradient(v, B, p0, max_iter, step_size, epsilon, plot_flag, p_opt_solver, fval_solver);
% disp(['Sub final gap: ', num2str(obj_values_sub(end))]);
disp(['Subgradient time: ', num2str(time_sub), ' seconds']);
disp(['Subgradient iterations: ', num2str(iter_sub)]);

%%% * - solve the problem by mirror descent
%%% ! Require changes for update and function value
[solution_md, time_md, iter_md, obj_values_md, distance_md] = quasi_primal_md(v, B, x0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver);
disp(['MD iterations: ', num2str(iter_md)]);
disp(['MD time: ', num2str(time_md), ' seconds']);

%%% * - solve the problem by single-agd
%%% ! Require changes for gradient and function value
plot_flag = true;
plot_flag_smooth = false;
adaptive = false;
delta_agd = 0.1; %%% Todo: AGD- smoothing parameter
[solution_agd, time_agd, iter_agd, obj_values_agd, dis_agd, convergence_agd] = quasi_dual_agd(v, B, mu0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta_agd, plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive);
disp(['AGD iterations: ', num2str(iter_agd)]);
disp(['AGD time: ', num2str(time_agd), ' seconds']);

%%% * - solve the problem by multiple-agd
%%% ! Require changes for gradient and function value
adaptive_plot_flag = true;  % Set to true if you want to plot the results
plot_flag = false;
plot_flag_smooth = false;
adaptive = true;
phase_num = 20; %%% Todo: SGR - parameter
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive] = quasi_dual_adaptive(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive, phase_num);
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);
disp(['Adaptive AGD time: ', num2str(total_time_adaptive), ' seconds']);

%%% * - plot the descent graph
x_md = 1:length(obj_values_md);
x_subgrad = 1:length(obj_values_sub);
x_agd = 1:length(obj_values_agd);
% x_adaptive = 1:length(obj_values_adaptive);

figure;
semilogy(x_md, abs(obj_values_md), '-d', 'DisplayName', 'Mirror Descent', 'LineWidth', 2); % Plot obj_mirror
hold on;
semilogy(x_subgrad, abs(obj_values_sub), '-d', 'DisplayName', 'Tatonnement', 'LineWidth', 2); % Plot obj_subgrad
semilogy(x_agd, abs(obj_values_agd), '-d', 'DisplayName', 'AGD', 'LineWidth', 2);
% semilogy(x_adaptive, abs(obj_values_adaptive), '-d', 'DisplayName', 'SGA', 'LineWidth', 2);
hold off;

% Set font sizes and other properties
set(gca, 'FontSize', 15); % Set axis font size
xlabel('Iteration', 'FontSize', 25); % X-axis label with larger font size
ylabel('Objective Value Gap', 'FontSize', 25); % Y-axis label with larger font size
title(['n=', num2str(n), ', m=', num2str(m)], 'FontSize', 25); % Graph title with larger font size
legend show; % Show legend
grid on; % Enable grid
