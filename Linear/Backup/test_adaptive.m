clc
clear
n = 10; % number of rows
m = 10; % number of columns
B = rand(n, 1); % Draw budgets from uniform distribution
B = B / sum(B); % Normalize to sum to 1
v = exprnd(1, n, m ); % Draw valuations from exponential distribution
v = v ./ sum(v, 2); % Normalize each row to sum to 1
plot_flag = true;

% Parameters for the algorithms
[x0, p0, mu_0, max_iter, step_size, eta, epsilon, L, sigma, mu_lower, mu_upper, delta] = linear_gen_par(v, B);

% Solve the problem by solver
p_opt_solver = ones(1,m); % Done by us - no actual use
 
% MD as solver
fprintf('Testing mirror descent...\n');
x_mirror = linear_primal_md(v, B, x0, eta, epsilon, max_iter,p_opt_solver, plot_flag);
p_mirror = sum(x_mirror, 1);
obj_mirror = sum(p_mirror) - sum(B .* log(min(p_mirror ./ v, [], 2)));
fprintf('Optimal Value by Solver', obj_mirror);

% Test the adaptive method for AGM
phase = 10;
epsilon_phase_change = 1e-1;
[solution, time, obj_values, f_smooth_values, dis_agd] = linear_dual_adaptive_agd(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, p_mirror, obj_mirror, phase, epsilon_phase_change);


%%% Todo: Compare algorithm - Totament, AGD and MD - Efficient AGD - Then Compare[solution, time, obj_values, f_smooth_values, dis_agd] = linear_dual_adaptive_agd(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, p_opt_solver, fval_solver, phase, epsilon_phase_change)