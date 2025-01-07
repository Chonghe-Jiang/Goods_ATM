clc
clear
%%% ! Add path
addpath('C:\Users\s1155203585\Dropbox\EG_EXP\Linear\EC_Submission\Efficiency');
n = 10; % number of rows
m = 10; % number of columns

% Generate budgets from uniform distribution and normalize
B = ones(n, 1); % Draw budgets from uniform distribution
B = B / sum(B); % Normalize to sum to 1

% Generate valuations from exponential distribution and normalize
v = exprnd(1, n, m); % Draw valuations from exponential distribution
v = v ./ sum(v, 2); % Normalize each row to sum to 1
plot_flag = true;
[p_opt, beta_opt, objective_value, solve_time] = linear_dual_solver(n, m, B, v);
% Compute the gap
mu = log(p_opt);
[gap, gap_array, matrix_backup] = linear_compute_gap(v, B, mu);
fprintf('Gap: %f\n', gap);

% % Parameters for the algorithms
% [x0, p0, mu_0, max_iter, step_size, eta, epsilon, L, sigma, mu_lower, mu_upper, delta] = linear_gen_par(v, B);
% p_opt_solver = ones(1, m); % Done by us - no actual use

% % MD as solver - To get the optimal point
% fprintf('Testing mirror descent...\n');
% [p_mirror, obj_mirror, x_mirror] = linear_primal_md(v, B, x0, eta, epsilon, max_iter, p_opt_solver, plot_flag);

