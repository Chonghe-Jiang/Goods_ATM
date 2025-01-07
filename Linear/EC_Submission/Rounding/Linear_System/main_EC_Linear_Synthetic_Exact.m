clc
clear
addpath('C:\Users\s1155203585\Dropbox\EG_EXP\Linear\EC_Submission\Efficiency');
%%% Todo: Consider the machine accuracy - important - whether the current stopping creteria in the inner problem is okay
%%% * Now the methodology is right and okay
n = 5;  % Number of rows
m = 5;   % Number of columns
B = ones(n, 1);  % Random B vector
v = rand(n,m);
v = v ./ sum(v, 2); 

%%%* Set common parameters for iterative methods
max_iter = 10000;
max_iter_adaptive = 4500;
p_lower = max(v .* B ./ sum(abs(v),2));
p_upper = norm(B, 1) * ones(1, m);
mu_lower = log(p_lower);
mu_upper = log(p_upper);
delta = 0.3;  
epsilon = 0.1;
sigma = min(exp(mu_lower));
L = exp(max(mu_upper)) + (sum(B) / delta); 

%%% * - ini of p0 and mu0
p0 = linear_init_gd(p_lower,p_upper,sum(B));
mu0 = log(p0); %
[p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
%%% *　Notice that the solver result is also not accurate 

%%% * To see the solution quality by the solver
% gap_solution = linear_compute_gap_cheating(v, B, log(p_opt_solver));
% disp(['Solution gap: ', num2str(gap_solution)]);
% disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Transfer to a row vector

adaptive_plot_flag = true;  % Set to true if you want to plot the results
plot_flag = false;
plot_flag_smooth = false;
adaptive = true;
phase_num = 20;
epsilon_current = 1e-4;
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive, results_matrix] = linear_dual_adaptive_exact(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive, phase_num, epsilon_current);
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);
disp(['Adaptive AGD time: ', num2str(total_time_adaptive), ' seconds']);

%%% * Extra test; See whether the optimal point satisfies thi
% mu_opt_solver = log(p_opt_solver);
% test_mat = log(v) - mu_opt_solver;
% max_values = max(test_mat, [], 2);
% binary_matrix = abs(test_mat - max_values) < 1e-4;
% test_result = linear_max_flow(mu_opt_solver,B,v,binary_matrix);
% disp('Test Result:');
% disp(test_result);

%%% * Output the final result
gap_adaptive = norm(solution_adaptive - log(p_opt_solver), 'fro');
disp(['Gap between Adaptive Solution and Optimal Solution: ', num2str(gap_adaptive)]);