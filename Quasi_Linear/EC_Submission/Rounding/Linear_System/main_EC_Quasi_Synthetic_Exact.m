%%% EC 2025 Submission version
clc
clear
addpath('/Users/chjiang/Dropbox/EG_EXP/Quasi_Linear/EC_Submission/Efficiency');

%%% ! Step 1: Basic setting of the problem
%%% Todo: Consider the machine accuracy 
n = 20;  % Number of rows
m = 20;   % Number of columns
B = ones(n, 1);  % Random B vector

% Define the folder name
dataset_folder = 'exact_data';

% Generate the filename for solver results based on n and m
solver_filename = sprintf('solver_exact_integer_%d_%d.mat', n, m);
solver_filepath = fullfile(dataset_folder, solver_filename); % Full path to the file
v_filename = sprintf('v_exact_integer_%d_%d.mat', n, m);
v_filepath = fullfile(dataset_folder, v_filename); % Full path to the file

% Check if the file exists. If it does, load 'v' from the file. Otherwise, generate 'v' and save it.
if exist(v_filepath, 'file') == 2
    load(v_filepath, 'v');  % Load 'v' from the file
    disp(['Loaded v from ', v_filepath]);
else
    v = randi([10,100],n, m); 
    % v = exprnd(10, n, m); does not work
    % v = lognrnd(0, 10, n, m); % Draw valuations from log-normal distribution
    % v = rand(n,m);
    % v = v ./ sum(v, 2); 
    save(v_filepath, 'v');  % Save 'v' to a file for future use
    disp(['Generated v and saved to ', v_filepath]);
end

max_iter = 10000;
max_iter_adaptive = 4500;
%%% ! Box constraint - Quasi Linear Version
p_lower = max(v .* B ./ (sum(abs(v),2)+B)); %%%!  quasi linear version
p_upper = max(v); %%%!  quasi linear version
mu_lower = log(p_lower);
mu_upper = log(p_upper);
delta = 0.3;  % Todo: parameter setting of iterative algorithm
epsilon = 0.1; %%% ! Note that this epsilon has no usage actually
sigma = min(exp(mu_lower));
L = exp(max(mu_upper)) + (sum(B) / delta); 
p0 = quasi_init_gd(p_lower, p_upper, sum(B));
mu0 = log(p0); %

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

%%% ! Step 2: Start the SGR-Exact: new version of this
adaptive_plot_flag = false;  %%% ! Set it to false version
plot_flag = false;
plot_flag_smooth = false;
adaptive = true;
phase_num = 20;
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive, results_matrix] = quasi_dual_adaptive_exact(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive, phase_num);
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);
disp(['Adaptive AGD time: ', num2str(total_time_adaptive), ' seconds']);
Exact_solver_gap = norm(solution_adaptive - log(p_opt_solver));
disp(['Solution Gap between Exact and Solver ', num2str(Exact_solver_gap)]);