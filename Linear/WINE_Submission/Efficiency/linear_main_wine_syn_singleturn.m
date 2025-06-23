%%% EC 2025 Submission
%%% ! Version: EC - Linear - Sythetic Data - Iterative method
%%% ! Optimal Version for Submission
%%% ! Includes Direct GD, NAG, MD, DA, Subgradient comparison with selectable algorithms
%%% ! Added integrated smoothed objective plot when plot_flag_smooth = true
%%% ! Adjusted plotting format and extended plots to max iteration count
clc
clear

% FOR DUAL AVERAGING, THE OPTIMAL PARAMETERS ARE:
% We only consider the exp generation
% 50*50 0.01
% 100*100 0.005
% 300*300 0.002

% The only thing we need to do is the multiturn experiment

% --- Algorithm Selection ---
run_subgradient = false;      % Set to true to run Subgradient (Tatonnement)
run_mirror_descent = true;   % Set to true to run Mirror Descent (Primal)
run_dual_averaging = false;    % Set to true to run Dual Averaging (Primal)
run_gd = false;              % Set to true to run Vanilla GD (Dual Smoothed)
run_agd = true;             % Set to true to run NAG (Dual Smoothed)
% run_projected_subgradient = false; % Removed
run_primal_pgd = false;        % Set to true to run Primal Projected Gradient Descent
run_primal_dual_pdhcg = false; % Set to true to run PDHCG
run_primal_dual_pdhg = true;   % Set to true to run PDHG (NEW)
% --------------------------

% Define problem parameters
n = 500;  % Number of rows
m = 500;   % Number of columns
B = ones(n,1); % For linear utilities, B_i = 1 is CEEI (w_i in PDHG paper)

% Define the folder name
dataset_folder = 'synthetica_dataset';

% Iteration and tolerance parameters
max_iter_sub = 20000; % Max iterations for Subgradient
max_iter_md = 20000;  % Max iterations for Mirror Descent
max_iter_da = 20000;  % Max iterations for Dual Averaging
max_iter_gd = 20000;  % Max iterations for Vanilla GD
max_iter_agd = 20000; % Max iterations for NAG (adjust as needed)
% max_iter_psub = 20000; % Removed
max_iter_primal_pgd = 20000; % Max iterations for Primal PGD

% PDHCG/PDHG specific parameters
max_outer_iter_pdhg = 100; % Number of outer loop restarts for PDHG/PDHCG
K_inner_iter_pdhg = 200; % Inner iterations per outer loop for PDHG/PDHCG
% Step sizes for PDHG/PDHCG (sigma and tau in paper). Need tuning.
% L is sqrt(n) for problem (9) in PDHCG paper (Lemma 1). For PDHG (problem 2) it's different.
% For PDHG (problem 2), L = sqrt(n) is for the (x,t) block to (p,y) block.
% The paper states L = sqrt(n).
sigma_pdhg = 5/(2*sqrt(n)); % Example initial value, requires tuning
tau_pdhg = 4/(2*sqrt(n)); % Example initial value, requires tuning


epsilon = 1e-1; % Stopping criteria based on original objective gap
delta = 0.1; % Smoothing parameter (used by GD and NAG)
eta_da = 0.5;  % Dual Averaging stepsize
plot_flag = false; % Set to true to see individual algorithm plots (Original Gap, Smooth Value, Distance)
plot_flag_smooth = true; % Set to true to see separate plot for smoothed objective value comparison (for GD and NAG)

% Generate the filename for solver results based on n and m
solver_filename = sprintf('solver_linear_rand_%d_%d.mat', n, m);
solver_filepath = fullfile(dataset_folder, solver_filename); % Full path to the file
v_filename = sprintf('v_linear_rand_%d_%d.mat', n, m);
v_filepath = fullfile(dataset_folder, v_filename); % Full path to the file

% Check if the file exists. If it does, load 'v' from the file. Otherwise, generate 'v' and save it.
if exist(v_filepath, 'file') == 2
    load(v_filepath, 'v');  % Load 'v' from the file
    disp(['Loaded v from ', v_filepath]);
else
    % v is u_ij in PDHG paper
    v = exprnd(10, n, m); % Generate valuations (Example: Exponential distribution)
    % Ensure v has no zero rows or columns for non-degeneracy (Assumption 1 in PDHG paper)
    if any(sum(v, 2) == 0)
        v(sum(v, 2) == 0, :) = exprnd(10, sum(sum(v, 2) == 0), m);
    end
    if any(sum(v, 1) == 0)
        v(:, sum(v, 1) == 0) = exprnd(10, n, sum(sum(v, 1) == 0));
    end
    v = v ./ sum(v, 2); % Normalize utilities (preprocessing step in PDHG paper, Sec 5.2)
    save(v_filepath, 'v');  % Save 'v' to a file for future use
    disp(['Generated v and saved to ', v_filepath]);
end

% Check if the solver results file exists. If it does, load the results. Otherwise, solve and save the results.
if exist(solver_filepath, 'file') == 2
    load(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Load solver results
    disp(['Loaded solver results from ', solver_filepath]);
else
    if exist('linear_dual_solver', 'file') == 2
        % fval_solver here is the optimal value of the *primal* problem (EG)
        % For linear utilities, the optimal primal and dual values are equal at equilibrium
        [p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
        save(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Save solver results
        disp(['Solved and saved solver results to ', solver_filepath]);
    else
        error('linear_dual_solver.m not found. Cannot compute optimal solution. Please ensure it is in your path.');
    end
end
disp(['Solver Optimal Value (Original Func): ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Ensure it's a row vector for consistency

%%% * Box constraint
% These bounds are relevant for dual methods that transform prices (log-prices)
p_lower = max(v .* B ./ sum(abs(v),2));
p_lower(p_lower <= 0) = 1e-12; % Replace non-positive with a small positive number
p_upper = norm(B, 1) * ones(1, m); % Upper bound for prices
mu_lower = log(p_lower);
mu_upper = log(p_upper);

%%% * Parameters
L_smooth = exp(max(mu_upper)) + (sum(B) / delta); % Lipschitz constant for smoothed gradient
step_size_sub = 1e-3; % Subgradient stepsize (needs tuning)
% step_size_psub = 1e-3; % Removed
step_size_primal_pgd = 1e-3; % Primal PGD stepsize
eta_md = 0.01;  % Mirror Descent stepsize
step_size_gd = 1 / L_smooth; % Step size for Vanilla GD

% Handle cases where sigma might be zero or negative
sigma_calc = min(exp(mu_lower));
if sigma_calc <=0
    warning('Calculated sigma is non-positive. Setting sigma to a small positive value.');
    sigma = 1e-10; % Assign a small positive value
else
    sigma = sigma_calc; % Strong convexity parameter (used by NAG)
end

adaptive = false; % Set to false for direct calls (enables/disables adaptive checks within functions)

%%% * - Initialization
if exist('linear_init_gd', 'file') ~= 2
    error('linear_init_gd.m not found. Cannot initialize dual variables.');
end
if exist('linear_init_md', 'file') ~= 2
    error('linear_init_md.m not found. Cannot initialize primal variables.');
end
% REMOVED: if exist('linear_init_primal_t', 'file') ~= 2
% REMOVED:     error('linear_init_primal_t.m not found. Cannot initialize auxiliary variable t.');
% REMOVED: end
% REMOVED: if exist('linear_init_dual_y', 'file') ~= 2
% REMOVED:     error('linear_init_dual_y.m not found. Cannot initialize dual variable y.');
% REMOVED: end


p0_dual_methods = linear_init_gd(p_lower,p_upper,sum(B)); % Initial feasible p for dual methods (p_0 for log-prices)
mu0 = log(p0_dual_methods); % Initial mu for dual methods (GD, NAG)
x0_primal_md = linear_init_md(p0_dual_methods,B); % Initial x for primal methods (MD, DA)

% Use x0_primal_md for Primal PGD and PDHG/PDHCG's x0
x0_for_primal_algs = x0_primal_md; 

% Use p0_dual_methods for PDHG/PDHCG's p0
p0_for_dual_space_algs = p0_dual_methods;

% *** MODIFICATION START: Initialize t0 and y0 directly in main function ***
t0_for_pdhg = ones(n, 1);
t0_for_pdhg(t0_for_pdhg <= 0) = 1e-12; % Ensure positivity t > 0

y0_for_pdhg = zeros(n, 1); % y can be initialized to zeros
% *** MODIFICATION END ***


disp('--- Running Selected Algorithms ---');

results = struct(); % Structure to store results of run algorithms
algo_iterations = []; % Store actual iterations for each run algorithm

%%% * - solve the problem by subgradient (Tatonnement) - Operates on Original Objective
if run_subgradient
    if exist('linear_dual_subgradient', 'file') == 2
        disp('Running Subgradient (Tatonnement)...');
        [solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = linear_dual_subgradient(v, B, p0_dual_methods, max_iter_sub, step_size_sub, epsilon, plot_flag, p_opt_solver, fval_solver);
        disp(['Subgradient time: ', num2str(time_sub), ' seconds, Iterations: ', num2str(iter_sub)]);
        results.subgradient.solution = solution_sub;
        results.subgradient.obj_values = obj_values_sub;
        results.subgradient.time = time_sub;
        results.subgradient.iter = iter_sub;
        results.subgradient.distance_final = norm(solution_sub - p_opt_solver);
        algo_iterations = [algo_iterations, iter_sub];
    else
        warning('linear_dual_subgradient.m not found. Skipping Subgradient.');
        run_subgradient = false;
    end
end

%%% * - solve the problem by Primal Projected Gradient Descent
if run_primal_pgd
    if exist('linear_primal_projected_gradient', 'file') == 2
        disp('Running Primal Projected Gradient Descent...');
        [solution_x_primal_pgd, solution_p_primal_pgd, obj_values_primal_pgd, dis_p_primal_pgd, time_primal_pgd, iter_primal_pgd] = linear_primal_projected_gradient(v, B, x0_for_primal_algs, max_iter_primal_pgd, step_size_primal_pgd, epsilon, plot_flag, p_opt_solver, fval_solver);
        disp(['Primal PGD time: ', num2str(time_primal_pgd), ' seconds, Iterations: ', num2str(iter_primal_pgd)]);
        results.primal_pgd.solution_x = solution_x_primal_pgd;
        results.primal_pgd.solution_p = solution_p_primal_pgd; % Store derived prices
        results.primal_pgd.obj_values = obj_values_primal_pgd;
        results.primal_pgd.time = time_primal_pgd;
        results.primal_pgd.iter = iter_primal_pgd;
        results.primal_pgd.distance_final = norm(solution_p_primal_pgd - p_opt_solver); % Compare prices
        algo_iterations = [algo_iterations, iter_primal_pgd];
    else
        warning('linear_primal_projected_gradient.m not found. Skipping Primal Projected Gradient Descent.');
        run_primal_pgd = false;
    end
end

%%% * - solve the problem by Primal-Dual PDHCG
if run_primal_dual_pdhcg
    if exist('linear_primal_dual_pdhcg', 'file') == 2
        disp('Running Primal-Dual PDHCG...');
        [solution_x_pdhcg, solution_p_pdhcg, obj_values_pdhcg, dis_p_pdhcg, time_pdhcg, iter_pdhcg] = linear_primal_dual_pdhcg(v, B, x0_for_primal_algs, p0_for_dual_space_algs, max_outer_iter_pdhg, K_inner_iter_pdhg, sigma_pdhg, tau_pdhg, epsilon, plot_flag, p_opt_solver, fval_solver);
        disp(['PDHCG time: ', num2str(time_pdhcg), ' seconds, Total Inner Iterations: ', num2str(iter_pdhcg)]);
        results.pdhcg.solution_x = solution_x_pdhcg;
        results.pdhcg.solution_p = solution_p_pdhcg; % Directly from algorithm
        results.pdhcg.obj_values = obj_values_pdhcg;
        results.pdhcg.time = time_pdhcg;
        results.pdhcg.iter = iter_pdhcg;
        results.pdhcg.distance_final = norm(solution_p_pdhcg - p_opt_solver); % Compare prices
        algo_iterations = [algo_iterations, iter_pdhcg];
    else
        warning('linear_primal_dual_pdhcg.m not found. Skipping Primal-Dual PDHCG.');
        run_primal_dual_pdhcg = false;
    end
end

%%% * - solve the problem by Primal-Dual PDHG (NEW)
if run_primal_dual_pdhg
    if exist('linear_primal_dual_pdhg', 'file') == 2
        disp('Running Primal-Dual PDHG...');
        % Pass x0, t0, p0, y0 directly
        [solution_x_pdhg, solution_t_pdhg, solution_p_pdhg, solution_y_pdhg, obj_values_pdhg, dis_p_pdhg, time_pdhg, iter_pdhg] = linear_primal_dual_pdhg(v, B, x0_for_primal_algs, t0_for_pdhg, p0_for_dual_space_algs, y0_for_pdhg, max_outer_iter_pdhg, K_inner_iter_pdhg, sigma_pdhg, tau_pdhg, epsilon, plot_flag, p_opt_solver, fval_solver);
        disp(['PDHG time: ', num2str(time_pdhg), ' seconds, Total Inner Iterations: ', num2str(iter_pdhg)]);
        results.pdhg.solution_x = solution_x_pdhg;
        results.pdhg.solution_t = solution_t_pdhg;
        results.pdhg.solution_p = solution_p_pdhg; % Directly from algorithm
        results.pdhg.solution_y = solution_y_pdhg;
        results.pdhg.obj_values = obj_values_pdhg;
        results.pdhg.time = time_pdhg;
        results.pdhg.iter = iter_pdhg;
        results.pdhg.distance_final = norm(solution_p_pdhg - p_opt_solver); % Compare prices
        algo_iterations = [algo_iterations, iter_pdhg];
    else
        warning('linear_primal_dual_pdhg.m not found. Skipping Primal-Dual PDHG.');
        run_primal_dual_pdhg = false;
    end
end


%%% * - solve the problem by mirror descent - Primal Method
if run_mirror_descent
     if exist('linear_primal_md', 'file') == 2
        disp('Running Mirror Descent...');
        [solution_md, time_md, iter_md, obj_values_md, distance_md] = linear_primal_md(v, B, x0_primal_md, eta_md, epsilon, max_iter_md, plot_flag, p_opt_solver, fval_solver);
        disp(['MD time: ', num2str(time_md), ' seconds, Iterations: ', num2str(iter_md)]);
        results.md.solution = solution_md;
        results.md.obj_values = obj_values_md;
        results.md.time = time_md;
        results.md.iter = iter_md;
        results.md.distance_final = norm(solution_md - p_opt_solver);
        algo_iterations = [algo_iterations, iter_md];
     else
        warning('linear_primal_md.m not found. Skipping Mirror Descent.');
        run_mirror_descent = false;
     end
end

%%% * - solve the problem by dual averaging - Primal Method
if run_dual_averaging
     if exist('linear_primal_da', 'file') == 2
        disp('Running Dual Averaging...');
        [solution_da, time_da, iter_da, obj_values_da, distance_da] = linear_primal_da(v, B, x0_primal_md, eta_da, epsilon, max_iter_da, plot_flag, p_opt_solver, fval_solver);
        disp(['DA time: ', num2str(time_da), ' seconds, Iterations: ', num2str(iter_da)]);
        results.da.solution = solution_da;
        results.da.obj_values = obj_values_da;
        results.da.time = time_da;
        results.da.iter = iter_da;
        results.da.distance_final = norm(solution_da - p_opt_solver);
        algo_iterations = [algo_iterations, iter_da];
     else
        warning('linear_primal_da.m not found. Skipping Dual Averaging.');
        run_dual_averaging = false;
     end
end

%%% * - solve the problem by Vanilla Gradient Descent (Direct Call) - Dual Smoothed
if run_gd
    if exist('linear_dual_gd', 'file') == 2
        disp('Running Vanilla GD...');
        [solution_gd, time_gd, iter_gd, obj_values_gd, f_smooth_values_gd, dis_gd, convergence_gd] = ...
            linear_dual_gd(v, B, mu0, max_iter_gd, L_smooth, epsilon, mu_lower, mu_upper, delta, ...
                           plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, step_size_gd, adaptive);
        disp(['Vanilla GD time: ', num2str(time_gd), ' seconds, Iterations: ', num2str(iter_gd)]);
        results.gd.solution = solution_gd;
        results.gd.obj_values = obj_values_gd;
        results.gd.f_smooth_values = f_smooth_values_gd;
        results.gd.time = time_gd;
        results.gd.iter = iter_gd;
        results.gd.distance_final = norm(exp(solution_gd) - p_opt_solver);
        algo_iterations = [algo_iterations, iter_gd];
    else
        warning('linear_dual_gd.m not found. Skipping Vanilla GD.');
        run_gd = false;
    end
end

%%% * - solve the problem by Nesterov's Accelerated Gradient (NAG - Direct Call) - Dual Smoothed
if run_agd
     if exist('linear_dual_agd', 'file') == 2
        disp('Running NAG...');
        [solution_agd, time_agd, iter_agd, obj_values_agd, f_smooth_values_agd, dis_agd, convergence_agd] = ...
            linear_dual_agd(v, B, mu0, max_iter_agd, L_smooth, sigma, epsilon, mu_lower, mu_upper, delta, ...
                            plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive);
        disp(['NAG time: ', num2str(time_agd), ' seconds, Iterations: ', num2str(iter_agd)]);
        results.agd.solution = solution_agd;
        results.agd.obj_values = obj_values_agd;
        results.agd.f_smooth_values = f_smooth_values_agd;
        results.agd.time = time_agd;
        results.agd.iter = iter_agd;
        results.agd.distance_final = norm(exp(solution_agd) - p_opt_solver);
        algo_iterations = [algo_iterations, iter_agd];
     else
        warning('linear_dual_agd.m not found or not modified to return smoothed values. Skipping NAG.');
        run_agd = false;
     end
end

disp('--- Algorithm Execution Finished ---');

% --- Find maximum common iteration count for plotting ---
max_common_iter = 0;
if ~isempty(algo_iterations)
    max_common_iter = max(algo_iterations);
    max_common_iter = max(1, max_common_iter);
    fprintf('Plotting results up to iteration: %d\n', max_common_iter);
else
    max_common_iter = 0;
end

%%% * - plot the combined descent graph (Using Original Objective Gap for comparison)
if max_common_iter > 0
    disp('Plotting Objective Gap results...');
    figure;
    hold on;
    legend_entries = {};
    plot_handles = [];
    all_obj_values_to_plot = [];

    plot_range_full = 1:max_common_iter;

    if run_mirror_descent && isfield(results, 'md') && ~isempty(results.md.obj_values)
        iters_actual = results.md.iter-1;
        plot_range_algo = 1:iters_actual;
        h = semilogy(plot_range_algo, abs(results.md.obj_values(plot_range_algo)), '-o', 'DisplayName', 'Mirror Descent (Primal Gap)', 'LineWidth', 2, 'MarkerSize', 4);
        plot_handles(end+1) = h;
        legend_entries{end+1} = 'Mirror Descent (Primal Gap)';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.md.obj_values(:))];
    end

    if run_dual_averaging && isfield(results, 'da') && ~isempty(results.da.obj_values)
        iters_actual = results.da.iter-1;
        plot_range_algo = 1:iters_actual;
        h = semilogy(plot_range_algo, abs(results.da.obj_values(plot_range_algo)), '-*', 'DisplayName', 'Dual Averaging (Primal Gap)', 'LineWidth', 2, 'MarkerSize', 4);
        plot_handles(end+1) = h;
        legend_entries{end+1} = 'Dual Averaging';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.da.obj_values(:))];
    end

    if run_subgradient && isfield(results, 'subgradient') && ~isempty(results.subgradient.obj_values)
        iters_actual = results.subgradient.iter-1; % Todo: tricky here
        plot_range_algo = 1:iters_actual;
        h = semilogy(plot_range_algo, abs(results.subgradient.obj_values(plot_range_algo)), '-s', 'DisplayName', 'Tatonnement (Dual Gap)', 'LineWidth', 2, 'MarkerSize', 4);
        plot_handles(end+1) = h;
        legend_entries{end+1} = 'Tatonnement ';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.subgradient.obj_values(:))];
    end
    
    if run_primal_pgd && isfield(results, 'primal_pgd') && ~isempty(results.primal_pgd.obj_values)
        iters_actual = results.primal_pgd.iter-1;
        plot_range_algo = 1:iters_actual;
        h = semilogy(plot_range_algo, abs(results.primal_pgd.obj_values(plot_range_algo)), '-^', 'DisplayName', 'Primal PGD (Primal Gap)', 'LineWidth', 2, 'MarkerSize', 4);
        plot_handles(end+1) = h;
        legend_entries{end+1} = 'Primal PGD ';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.primal_pgd.obj_values(:))];
    end

    if run_primal_dual_pdhcg && isfield(results, 'pdhcg') && ~isempty(results.pdhcg.obj_values)
        iters_actual = results.pdhcg.iter-1;
        plot_range_algo = 1:iters_actual;
        h = semilogy(plot_range_algo, abs(results.pdhcg.obj_values(plot_range_algo)), '-v', 'DisplayName', 'Primal-Dual PDHCG', 'LineWidth', 2, 'MarkerSize', 4);
        plot_handles(end+1) = h;
        legend_entries{end+1} = 'Primal-Dual PDHCG';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.pdhcg.obj_values(:))];
    end

    if run_primal_dual_pdhg && isfield(results, 'pdhg') && ~isempty(results.pdhg.obj_values)
        iters_actual = results.pdhg.iter-1;
        plot_range_algo = 1:iters_actual;
        h = semilogy(plot_range_algo, abs(results.pdhg.obj_values(plot_range_algo)), '-<', 'DisplayName', 'Primal-Dual PDHG', 'LineWidth', 2, 'MarkerSize', 4);
        plot_handles(end+1) = h;
        legend_entries{end+1} = 'Primal-Dual PDHG';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.pdhg.obj_values(:))];
    end


    if run_gd && isfield(results, 'gd') && ~isempty(results.gd.obj_values)
        iters_actual = results.gd.iter-1;
        plot_range_algo = 1:iters_actual;
        h = semilogy(plot_range_algo, abs(results.gd.obj_values(plot_range_algo)), '-.', 'DisplayName', 'Vanilla GD (Dual Gap)', 'LineWidth', 2, 'MarkerSize', 4);
        plot_handles(end+1) = h;
        legend_entries{end+1} = 'Vanilla GD ';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.gd.obj_values(:))];
    end

    if run_agd && isfield(results, 'agd') && ~isempty(results.agd.obj_values)
        iters_actual = results.agd.iter-1;
        plot_range_algo = 1:iters_actual;
        h = semilogy(plot_range_algo, abs(results.agd.obj_values(plot_range_algo)), '-d', 'DisplayName', 'NAG (Dual Gap)', 'LineWidth', 2, 'MarkerSize', 4);
        plot_handles(end+1) = h;
        legend_entries{end+1} = 'NAG';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.agd.obj_values(:))];
    end

    hold off;

    set(gca, 'FontSize', 15);
    xlabel('Iteration', 'FontSize', 20);
    ylabel('Objective Value Gap F(\cdot) - F*', 'FontSize', 20);
    title_str_obj = sprintf('Objective Gap Comparison (n=%d, m=%d)', n, m);
    title(title_str_obj, 'FontSize', 25);
    if ~isempty(legend_entries)
        lgd = legend(plot_handles, legend_entries, 'Location', 'best');
        lgd.FontSize = 15;
    else
        warning('No algorithms produced results for objective gap plotting.');
    end
    grid on;
    set(gca, 'YScale', 'log');
    xlim([1, max_common_iter]);

    valid_obj_values = all_obj_values_to_plot(all_obj_values_to_plot > 0 & isfinite(all_obj_values_to_plot));
    if ~isempty(valid_obj_values)
        min_val = max(min(valid_obj_values)*0.1, epsilon/10);
        max_val = max(valid_obj_values)*10;
        if min_val <= 0 || ~isfinite(min_val)
            min_val = epsilon/100;
        end
        if max_val <= min_val || ~isfinite(max_val)
            max_val = max(1, min_val * 100);
        end
        ylim([min_val, max_val]);
    else
        ylim([epsilon/10, 1]);
    end
else
    disp('No algorithms selected or ran successfully for objective gap plot.');
end

%%% * -- Plot the integrated smoothed objective function graph -- %%%
% Todo: we do not need smoothness here anymore
% if plot_flag_smooth && run_gd && run_agd && max_common_iter > 0
%     if isfield(results, 'gd') && isfield(results.gd, 'f_smooth_values') && ~isempty(results.gd.f_smooth_values) && ...
%        isfield(results, 'agd') && isfield(results.agd, 'f_smooth_values') && ~isempty(results.agd.f_smooth_values)

%         disp('Plotting Smoothed Objective Function Comparison...');
%         figure;
%         hold on;
%         plot_handles_smooth = [];

%         if run_gd && isfield(results, 'gd')
%             iters_gd = results.gd.iter;
%             plot_range_gd = 1:iters_gd;
%             h = plot(plot_range_algo, results.gd.f_smooth_values(plot_range_algo), '-.', 'DisplayName', 'Vanilla GD', 'LineWidth', 2, 'MarkerSize', 4);
%             plot_handles_smooth(end+1) = h;
%         end

%         if run_agd && isfield(results, 'agd')
%             iters_agd = results.agd.iter;
%             plot_range_agd = 1:iters_agd;
%             h = plot(plot_range_algo, results.agd.f_smooth_values(plot_range_algo), '-d', 'DisplayName', 'NAG', 'LineWidth', 2, 'MarkerSize', 4);
%             plot_handles_smooth(end+1) = h;
%         end

%         hold off;

%         set(gca, 'FontSize', 15);
%         xlabel('Iteration', 'FontSize', 25);
%         ylabel('Smoothed Objective Value F_{\delta}(\mu)', 'FontSize', 25);
%         title_str_smooth = sprintf('Comparison on Problem (P_{\\delta}) (n=%d, m=%d)', n, m);
%         title(title_str_smooth, 'FontSize', 25);
%         lgd_smooth = legend(plot_handles_smooth, 'Location', 'best');
%         lgd_smooth.FontSize = 15;
%         grid on;
%         xlim([1, max_common_iter]);

%     else
%         warning('plot_flag_smooth is true, but smoothed objective results for both GD and NAG are not available. Skipping plot.');
%     end
% elseif plot_flag_smooth
%     disp('plot_flag_smooth is true, but both GD and NAG were not selected or did not run successfully, or no common iterations found. Skipping smoothed objective plot.');
% end


%%% * - Calculate and Print final distances to optimal primal solution p*
disp('--- Final Distance to Optimal Primal Solution p* ---');
found_results = false;

if run_mirror_descent && isfield(results, 'md')
    disp(['Mirror Descent: ', num2str(results.md.distance_final)]);
    found_results = true;
end

if run_dual_averaging && isfield(results, 'da')
    disp(['Dual Averaging: ', num2str(results.da.distance_final)]);
    found_results = true;
end

if run_subgradient && isfield(results, 'subgradient')
    disp(['Tatonnement (Subgradient): ', num2str(results.subgradient.distance_final)]);
    found_results = true;
end

if run_primal_pgd && isfield(results, 'primal_pgd')
    disp(['Primal PGD: ', num2str(results.primal_pgd.distance_final)]);
    found_results = true;
end

if run_primal_dual_pdhcg && isfield(results, 'pdhcg')
    disp(['Primal-Dual PDHCG: ', num2str(results.pdhcg.distance_final)]);
    found_results = true;
end

if run_primal_dual_pdhg && isfield(results, 'pdhg')
    disp(['Primal-Dual PDHG: ', num2str(results.pdhg.distance_final)]);
    found_results = true;
end

if run_gd && isfield(results, 'gd')
    disp(['Vanilla GD: ', num2str(results.gd.distance_final)]);
    found_results = true;
end

if run_agd && isfield(results, 'agd')
    disp(['NAG: ', num2str(results.agd.distance_final)]);
    found_results = true;
end

if ~found_results
    disp('No results to display.');
end

disp('--- Script Finished ---');