%%% EC 2025 Submission
%%% ! Version: EC - Quasi Linear - Synthetic Data - Iterative method
%%% ! Optimal Version for Submission
%%% ! In QL Setting -> we do not do normalization anymore
%%% ! REVISED PLOTTING with user-specified RGB colors
clc
clear

% Define problem parameters
n = 500;  % Number of rows 
m = 500;   % Number of columns
B = ones(n, 1);  % Random B vector

% --- Iteration and Tolerance Parameters ---
epsilon = 1e-4; 
max_iter = 20000;
max_iter_adaptive_agd = 4000; % Max iterations per phase for APM
phase_num = 20; %%% Todo: SGR - parameter for APM
switch_step_sub = 10000; % Iteration count to switch step size in Adaptive Tatonnement
switch_step_md = 500; % Iteration count to switch step size in Adaptive Mirror Descent

% --- Algorithm Selection ---
% In this script, we run all algorithms for comparison
run_subgradient = true;      % Tatonnement
run_adaptive_sub = true;   % Adaptive Tatonnement
run_primal_md = true;      % Mirror Descent
run_adaptive_md = true;    % Adaptive Mirror Descent
run_adaptive_nag = true;     % APM (Adaptive Accelerated Proximal Method)

% --- Plotting Flags ---
plot_overall_adaptive = false; % Master flag to enable/disable plotting within functions

% Define the folder name
dataset_folder = 'synthetica_dataset_quasi';  % Changed folder name to include "quasi"
% Generate the filename for v based on n and m
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
    % v = exprnd(10, n, m);
    v = exprnd(10, n, m);  % Generate a random matrix of size n x m 
    % v = v ./ sum(v, 2); 
    save(v_filepath, 'v');  % Save 'v' to a file for future use
    disp(['Generated v and saved to ', v_filepath]);
end

% Check if the solver results file exists. If it does, load the results. Otherwise, solve and save the results.
if exist(solver_filepath, 'file') == 2
    load(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Load solver results from the file
    disp(['Loaded solver results from ', solver_filepath]);
else
    % Solve the problem using the solver
    if exist('quasi_dual_solver', 'file') == 2
        [p_opt_solver, beta_opt, fval_solver, solve_time] = quasi_dual_solver(n, m, B, v);
        save(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Save solver results to a file
        disp(['Solved and saved solver results to ', solver_filepath]);
    else
        error('quasi_dual_solver.m not found. Cannot compute optimal solution.');
    end
end
disp(['Solver Optimal Value (Original Func): ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Transfer to a row vector

%%% ! Box constraint - Quasi Linear Version
p_lower = max(v .* B ./ (sum(abs(v), 2) + B)); %%%! quasi linear version
p_upper = max(v); %%%! quasi linear version
mu_lower = log(p_lower);
mu_upper = log(p_upper);

%%% * Parameter for convexity and smoothness and stepsize
delta = 0.1;  
sigma = min(exp(mu_lower));
L = exp(max(mu_upper)) + (sum(B) / delta); 
step_size_sub = 1e-5; %%% Todo: Subgradient stepsize
eta_md = 0.2;  %%% Todo: MD stepsize - now is the optimal step size 1

%%% * - Initialize p0 and mu0
if exist('quasi_init_gd', 'file') ~= 2 || exist('quasi_init_md', 'file') ~= 2
    error('Initialization functions (quasi_init_gd.m or quasi_init_md.m) not found.');
end
p0 = quasi_init_gd(p_lower, p_upper, sum(B));
mu0 = log(p0); 
x0 = quasi_init_md(p0, B);

disp('--- Running Selected Algorithms ---');

results = struct(); % Structure to store results

%%% * - Solve the problem by Tatonnement (Subgradient)
if run_subgradient
    if exist('quasi_dual_subgradient', 'file') == 2
        disp('Running Tatonnement...');
        [solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = ...
            quasi_dual_subgradient(v, B, p0, max_iter, step_size_sub, epsilon, plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['Tatonnement time: ', num2str(time_sub), ' seconds, Iterations: ', num2str(iter_sub)]);
        results.subgradient.solution = solution_sub;
        results.subgradient.obj_values = obj_values_sub;
        results.subgradient.time = time_sub;
        results.subgradient.iter = iter_sub;
        results.subgradient.distance_final = norm(solution_sub - p_opt_solver);
    else
        warning('quasi_dual_subgradient.m not found. Skipping Tatonnement.');
        run_subgradient = false;
    end
end

%%% * - Solve the problem by Adaptive Tatonnement (Adaptive Subgradient)
if run_adaptive_sub
    if exist('quasi_dual_subgradient_adaptive', 'file') == 2
        disp('Running Adaptive Tatonnement...');
        [solution_sub_adaptive, obj_values_sub_adaptive, dis_sub_adaptive, time_sub_adaptive, iter_sub_adaptive] = ...
            quasi_dual_subgradient_adaptive(v, B, p0, max_iter, step_size_sub, epsilon, plot_overall_adaptive, p_opt_solver, fval_solver, switch_step_sub);
        disp(['Adaptive Tatonnement time: ', num2str(time_sub_adaptive), ' seconds, Iterations: ', num2str(iter_sub_adaptive)]);
        results.adaptive_sub.solution = solution_sub_adaptive;
        results.adaptive_sub.obj_values = obj_values_sub_adaptive;
        results.adaptive_sub.time = time_sub_adaptive;
        results.adaptive_sub.iter = iter_sub_adaptive;
        results.adaptive_sub.distance_final = norm(solution_sub_adaptive - p_opt_solver);
    else
        warning('quasi_dual_subgradient_adaptive.m not found. Skipping Adaptive Tatonnement.');
        run_adaptive_sub = false;
    end
end

%%% * - Solve the problem by Mirror Descent
if run_primal_md
    if exist('quasi_primal_md', 'file') == 2
        disp('Running Mirror Descent...');
        [solution_md, time_md, iter_md, obj_values_md, distance_md] = ...
            quasi_primal_md(v, B, x0, eta_md, epsilon, max_iter, plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['Mirror Descent time: ', num2str(time_md), ' seconds, Iterations: ', num2str(iter_md)]);
        results.primal_md.solution = solution_md;
        results.primal_md.obj_values = obj_values_md;
        results.primal_md.time = time_md;
        results.primal_md.iter = iter_md;
        results.primal_md.distance_final = norm(solution_md - p_opt_solver);
    else
        warning('quasi_primal_md.m not found. Skipping Mirror Descent.');
        run_primal_md = false;
    end
end

%%% * - Solve the problem by Adaptive Mirror Descent
if run_adaptive_md
    if exist('quasi_primal_md_adaptive', 'file') == 2
        disp('Running Adaptive Mirror Descent...');
        [solution_md_adaptive, time_md_adaptive, iter_md_adaptive, obj_values_md_adaptive, distance_md_adaptive] = ...
            quasi_primal_md_adaptive(v, B, x0, epsilon, max_iter, plot_overall_adaptive, p_opt_solver, fval_solver, switch_step_md);
        disp(['Adaptive MD time: ', num2str(time_md_adaptive), ' seconds, Iterations: ', num2str(iter_md_adaptive)]);
        results.adaptive_md.solution = solution_md_adaptive;
        results.adaptive_md.obj_values = obj_values_md_adaptive;
        results.adaptive_md.time = time_md_adaptive;
        results.adaptive_md.iter = iter_md_adaptive;
        results.adaptive_md.distance_final = norm(solution_md_adaptive - p_opt_solver);
    else
        warning('quasi_primal_md_adaptive.m not found. Skipping Adaptive Mirror Descent.');
        run_adaptive_md = false;
    end
end

%%% * - Solve the problem by APM (Adaptive AGD)
if run_adaptive_nag
    if exist('quasi_dual_adaptive', 'file') == 2
        disp('Running APM...');
        adaptive_plot_flag = true; 
        plot_flag_internal = false;
        plot_flag_smooth = false;
        adaptive_flag_internal = true;
        [solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive] = ...
            quasi_dual_adaptive(v, B, mu0, max_iter_adaptive_agd, L, sigma, epsilon, mu_lower, mu_upper, delta, ...
                                plot_flag_internal, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive_flag_internal, phase_num);
        disp(['APM time: ', num2str(total_time_adaptive), ' seconds, Total Iterations: ', num2str(total_iter_adaptive)]);
        results.adaptive_nag.solution = solution_adaptive;
        results.adaptive_nag.obj_values = obj_values_adaptive;
        results.adaptive_nag.time = total_time_adaptive;
        results.adaptive_nag.iter = total_iter_adaptive;
        results.adaptive_nag.distance_final = norm(exp(solution_adaptive) - p_opt_solver);
    else
        warning('quasi_dual_adaptive.m not found. Skipping APM.');
        run_adaptive_nag = false;
    end
end

disp('--- Algorithm Execution Finished ---');

%%% * - plot the combined descent graph (Using Original Objective Gap for comparison)
if run_adaptive_nag || run_adaptive_sub || run_adaptive_md || run_primal_md || run_subgradient
    disp('Plotting Objective Gap results...');
    figure; % Create figure for comparison plot
    hold on;
    legend_entries = {}; % Cell array to store legend entries dynamically
    all_obj_values_to_plot = []; % Collect objective values for ylim calculation

    % Get the maximum iteration count across all run algorithms for consistent x-axis
    max_overall_iters = 0;
    if run_adaptive_nag && isfield(results, 'adaptive_nag'), max_overall_iters = max(max_overall_iters, results.adaptive_nag.iter); end
    if run_adaptive_sub && isfield(results, 'adaptive_sub'), max_overall_iters = max(max_overall_iters, results.adaptive_sub.iter); end
    if run_adaptive_md && isfield(results, 'adaptive_md'), max_overall_iters = max(max_overall_iters, results.adaptive_md.iter); end
    if run_primal_md && isfield(results, 'primal_md'), max_overall_iters = max(max_overall_iters, results.primal_md.iter); end
    if run_subgradient && isfield(results, 'subgradient'), max_overall_iters = max(max_overall_iters, results.subgradient.iter); end

    % Plot Tatonnement results
    if run_subgradient && isfield(results, 'subgradient') && ~isempty(results.subgradient.obj_values)
        total_iters = results.subgradient.iter;
        plot_range_algo = 1:(total_iters);
        if ~isempty(plot_range_algo) && total_iters > 0
            semilogy(plot_range_algo, abs(results.subgradient.obj_values(plot_range_algo)), '-v', 'DisplayName', 'Tatonnement', 'LineWidth', 2, 'MarkerSize', 4, 'Color', [0.00, 0.45, 0.74]);
            legend_entries{end+1} = 'Tatonnement';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.subgradient.obj_values(:))];
        end
    end
    
    % Plot Adaptive Tatonnement results
    if run_adaptive_sub && isfield(results, 'adaptive_sub') && ~isempty(results.adaptive_sub.obj_values)
        total_iters = results.adaptive_sub.iter;
        plot_range_algo = 1:(total_iters);
        if ~isempty(plot_range_algo) && total_iters > 0
            semilogy(plot_range_algo, abs(results.adaptive_sub.obj_values(plot_range_algo)), '-s', 'DisplayName', 'Adaptive Tatonnement', 'LineWidth', 2, 'MarkerSize', 4, 'Color', [0.85, 0.33, 0.10]);
            legend_entries{end+1} = 'Adaptive Tatonnement';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_sub.obj_values(:))];
        end
    end
    
    % Plot Mirror Descent results
    if run_primal_md && isfield(results, 'primal_md') && ~isempty(results.primal_md.obj_values)
        total_iters = results.primal_md.iter;
        plot_range_algo = 1:(total_iters-1);
        if ~isempty(plot_range_algo) && total_iters > 0
            semilogy(plot_range_algo, abs(results.primal_md.obj_values(plot_range_algo)), '-x', 'DisplayName', 'Mirror Descent', 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0.93, 0.69, 0.13]);
            legend_entries{end+1} = 'Mirror Descent';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.primal_md.obj_values(:))];
        end
    end

    % Plot Adaptive MD results
    if run_adaptive_md && isfield(results, 'adaptive_md') && ~isempty(results.adaptive_md.obj_values)
        total_iters = results.adaptive_md.iter;
        plot_range_algo = 1:(total_iters-1);
        if ~isempty(plot_range_algo) && total_iters > 0
            semilogy(plot_range_algo, abs(results.adaptive_md.obj_values(plot_range_algo)), '-o', 'DisplayName', 'Adaptive Mirror Descent', 'LineWidth', 2, 'MarkerSize', 4, 'Color', [0.49, 0.18, 0.56]);
            legend_entries{end+1} = 'Adaptive MD';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_md.obj_values(:))];
        end
    end

    % Plot APM results
    if run_adaptive_nag && isfield(results, 'adaptive_nag') && ~isempty(results.adaptive_nag.obj_values)
        total_iters = results.adaptive_nag.iter;
        plot_range_algo = 1:(total_iters - 1);
        if ~isempty(plot_range_algo) && total_iters > 0
            semilogy(plot_range_algo, abs(results.adaptive_nag.obj_values(plot_range_algo)), '-d', 'DisplayName', 'APM', 'LineWidth', 2, 'MarkerSize', 4, 'Color', [0.30, 0.75, 0.93]);
            legend_entries{end+1} = 'APM';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_nag.obj_values(:))];
        end
    end

    hold off;

    % Customize plot
    set(gca, 'FontSize', 15); % Axis font size
    xlabel('Iteration', 'FontSize', 20);
    ylabel('Objective Value Gap', 'FontSize', 20);
    title_str_obj = sprintf('n=%d, m=%d', n, m);
    title(title_str_obj, 'FontSize', 25);
    
    % Set logarithmic y-axis with proper ticks
    set(gca, 'YScale', 'log');
    grid on;
    
    % Dynamic y-axis limits and ticks
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
        
        min_exp = floor(log10(min_val));
        max_exp = ceil(log10(max_val));
        yticks(10.^(min_exp:1:max_exp));
        yticklabels(arrayfun(@(x) sprintf('10^{%d}', x), min_exp:1:max_exp, 'UniformOutput', false));
    else
        ylim([epsilon/10, 1]);
        yticks(10.^(-2:1:0));
        yticklabels({'10^{-2}', '10^{-1}', '10^{0}'});
    end
    
    % Set x-axis limits
    max_plot_iters = max_overall_iters;
    if max_plot_iters > 0
        xlim([1, max_plot_iters]);
    else
        xlim([1, 100]);
    end

    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'northeast', 'FontSize', 10);
    else
        warning('No algorithms produced results for objective gap plotting.');
    end
else
    disp('No algorithms selected or ran successfully for objective gap plot.');
end

% --- Final Distance to Optimal Primal Solution p* ---
disp('--- Final Distance to Optimal Primal Solution p* ---');
found_results = false;

if run_subgradient && isfield(results, 'subgradient')
    disp(['Tatonnement: ', num2str(results.subgradient.distance_final)]);
    found_results = true;
end

if run_adaptive_sub && isfield(results, 'adaptive_sub')
    disp(['Adaptive Tatonnement: ', num2str(results.adaptive_sub.distance_final)]);
    found_results = true;
end

if run_primal_md && isfield(results, 'primal_md')
    disp(['Mirror Descent: ', num2str(results.primal_md.distance_final)]);
    found_results = true;
end

if run_adaptive_md && isfield(results, 'adaptive_md')
    disp(['Adaptive Mirror Descent: ', num2str(results.adaptive_md.distance_final)]);
    found_results = true;
end

if run_adaptive_nag && isfield(results, 'adaptive_nag')
    disp(['APM: ', num2str(results.adaptive_nag.distance_final)]);
    found_results = true;
end

if ~found_results
    disp('No results to display.');
end

disp('--- Script Finished ---');