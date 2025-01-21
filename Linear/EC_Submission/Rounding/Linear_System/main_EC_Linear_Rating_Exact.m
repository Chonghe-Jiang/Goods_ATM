%%% EC 2025 Submission
%%% ! General
clc
clear
addpath('/Users/chjiang/Dropbox/EG_EXP/Linear/EC_Submission/Efficiency');
v = readmatrix('/Users/chjiang/Dropbox/EG_EXP/Linear/EC_Submission/Efficiency/Dataset/Ratings_kroer.csv');
v = v(1:100,1:100);
% v = v/norm(v,"fro");
%%% ! Step -1: whether to slicing the matrix or not
v = floor(v)+1;
[n,m] = size(v);
B = ones(n,1);
max_iter = 10000;
max_iter_adaptive = 4500;
p_lower = max(v .* B ./ sum(abs(v),2));
p_upper = norm(B, 1) * ones(1, m);
mu_lower = log(p_lower);
mu_upper = log(p_upper);
delta = 0.1;  % Todo: parameter setting of iterative algorithm
epsilon = 0.1; %%% ! Note that this epsilon has no usage actually
sigma = min(exp(mu_lower));
L = exp(max(mu_upper)) + (sum(B) / delta); 
p0 = linear_init_gd(p_lower, p_upper, sum(B));
mu0 = log(p0); %
%%% ! Step 0: Decide whether solve or not
[p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
% disp(['Solver optimal value: ', num2str(fval_solver)]);
% disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Transfer to a row vector

%%% ! Step 1: Load the Solver Results
% load('/Users/chjiang/Dropbox/EG_EXP/Linear/EC_Submission/Efficiency/Solver_Rating.mat', 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');
%%% ! Step 2: Calculate the Gap
gap_optimal = linear_optimal_gap(v, p_opt_solver);
disp(['The Gap is ', num2str(gap_optimal)]);
%%% ! Step 3: Start the SGR-Exact: new version of this
adaptive_plot_flag = false;  %%% ! Set it to false version
plot_flag = false;
plot_flag_smooth = false;
adaptive = true;
phase_num = 1000;
tic;
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive, results_matrix] = linear_dual_adaptive_exact(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive, phase_num);
adaptive_time = toc;
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);
disp(['Adaptive AGD time: ', num2str(adaptive_time), ' seconds']);
Exact_solver_gap = norm(solution_adaptive - log(p_opt_solver));
disp(['Solution Gap between Exact and Solver ', num2str(Exact_solver_gap)]);




