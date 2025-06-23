%%% EC 2025 Submission version
clc
clear
addpath('/Users/chjiang/Dropbox/EG_EXP/Linear/EC_Submission/Efficiency');

%%% ! Step 1: Basic setting of the problem
%%% Todo: Consider the machine accuracy 
n = 50;  % Number of rows
m = 50;   % Number of columns
B = ones(n, 1);  % Random B vector

%%% Todo: Change the data type here
% Define the folder name
dataset_folder = 'exact_data_EC';
% Define the CSV filename
csv_filename = 'solver_adaptive_gap_integer.csv';

% Generate the filename for solver results based on n and m
%%% Todo: Change the file name here
solver_filename = sprintf('solver_exact_integer_%d_%d.mat', n, m);
solver_filepath = fullfile(dataset_folder, solver_filename); % Full path to the file
v_filename = sprintf('v_exact_integer_%d_%d.mat', n, m);
v_filepath = fullfile(dataset_folder, v_filename); % Full path to the file

% Check if the file exists. If it does, load 'v' from the file. Otherwise, generate 'v' and save it.
if exist(v_filepath, 'file') == 2
    load(v_filepath, 'v');  % Load 'v' from the file
    disp(['Loaded v from ', v_filepath]);
else
    %%% Todo: Change the data type here
    % v = 10*rand(n, m); 
    % v = exprnd(10, n, m); 
    % v = lognrnd(0, 1, n, m);
    v = randi([1,10],n,m);
    save(v_filepath, 'v');  % Save 'v' to a file for future use
    disp(['Generated v and saved to ', v_filepath]);
end

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

% Check if the solver results file exists. If it does, load the results. Otherwise, solve and save the results.
if exist(solver_filepath, 'file') == 2
    load(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Load solver results from the file
    disp(['Loaded solver results from ', solver_filepath]);
else
    % Solve the problem using the solver
    [p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
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
phase_num = 1000;
tic;
[solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive, results_matrix] = linear_dual_adaptive_exact(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive, phase_num);
adaptive_time = toc;
disp(['Adaptive AGD iterations: ', num2str(total_iter_adaptive)]);
disp(['Adaptive AGD time: ', num2str(adaptive_time), ' seconds']);
Exact_solver_gap = norm(solution_adaptive - log(p_opt_solver));
disp(['Solution Gap between Exact and Solver ', num2str(Exact_solver_gap)]);
gap_optimal = linear_optimal_gap(v, solution_adaptive);
disp(['The Gap is ', num2str(gap_optimal)]);
%%% ! Step 3: Document n, m, solver_time, adaptive_time, and gap_optimal to a CSV file
% Create a table with the data
data_table = table(n, m, solve_time, adaptive_time, gap_optimal, ...
    'VariableNames', {'n', 'm', 'solver_time', 'adaptive_time', 'gap_optimal'});



% Check if the file exists. If it does, append the data. Otherwise, create a new file.
if exist(csv_filename, 'file') == 2
    % File exists, append the data
    writetable(data_table, csv_filename, 'WriteMode', 'append');
else
    % File does not exist, create a new file with headers
    writetable(data_table, csv_filename);
end

disp(['Data saved to ', csv_filename]);