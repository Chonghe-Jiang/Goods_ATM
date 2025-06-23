%%% EC 2025 Submission
%%% ! Version: EC - Linear - Sythetic Data - Multi-Turn/Adaptive Iterative methods
%%% ! Selectable Adaptive Algorithms Comparison
clc
clear

% --- Algorithm Selection ---
run_adaptive_gd = false;   % Set to true to run Adaptive Vanilla GD
run_adaptive_nag = true;  % Set to true to run Adaptive NAG (using linear_dual_adaptive)
run_adaptive_sub = false; % Set to true to run Adaptive Subgradient
run_adaptive_md = false;  % Set to true to run Adaptive Mirror Descent
run_dual_averaging = false; % Set to true to run Dual Averaging (non-adaptive version)
% --------------------------

% Define problem parameters
n = 50;  % Number of rows
m = 50;   % Number of columns
B = ones(n,1);

% Define the folder name
dataset_folder = 'synthetica_dataset';

% --- Iteration and Tolerance Parameters ---
% General
epsilon = 1e-4; % Stopping criteria based on original objective gap

% Adaptive GD / NAG Parameters
max_iter_per_phase = 1000; % Max iterations *per phase* for Adaptive GD/NAG
phase_num = 20;          % Number of phases for Adaptive GD/NAG
delta_initial = 0.1;     % Initial smoothing parameter (used by Adaptive GD/NAG)

% Adaptive Subgradient / MD Parameters
max_iter_adaptive_sub = 20000; % Max total iterations for Adaptive Subgradient
max_iter_adaptive_md = 20000;  % Max total iterations for Adaptive MD
switch_step_sub = 5000;       % Iteration count to switch step size in Adaptive Subgradient
switch_step_md = 5000;        % Iteration count to switch step size in Adaptive MD

% Dual Averaging Parameters
max_iter_da = 20000;          % Max iterations for Dual Averaging

% Plotting Flags (Control plotting *within* the called functions, if implemented)
plot_individual_phases = false; % e.g., plot for each GD/NAG phase call
plot_overall_adaptive = false; % e.g., plot overall convergence from within adaptive functions

% --- End Parameters ---

% Generate the filename for solver results based on n and m
solver_filename = sprintf('solver_linear_exp_%d_%d.mat', n, m);
solver_filepath = fullfile(dataset_folder, solver_filename); % Full path to the file
v_filename = sprintf('v_linear_exp_%d_%d.mat', n, m);
v_filepath = fullfile(dataset_folder, v_filename); % Full path to the file

% Check if the file exists. If it does, load 'v' from the file. Otherwise, generate 'v' and save it.
if exist(v_filepath, 'file') == 2
    load(v_filepath, 'v');  % Load 'v' from the file
    disp(['Loaded v from ', v_filepath]);
else
    v = exprnd(10, n, m); % Generate valuations
    v = v ./ sum(v, 2);
    save(v_filepath, 'v');  % Save 'v' to a file for future use
    disp(['Generated v and saved to ', v_filepath]);
end

% Check if the solver results file exists. If it does, load the results. Otherwise, solve and save the results.
if exist(solver_filepath, 'file') == 2
    load(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Load solver results
    disp(['Loaded solver results from ', solver_filepath]);
else
    if exist('linear_dual_solver', 'file') == 2
        [p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
        save(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');
        disp(['Solved and saved solver results to ', solver_filepath]);
    else
        error('linear_dual_solver.m not found. Cannot compute optimal solution.');
    end
end
disp(['Solver Optimal Value (Original Func): ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Ensure it's a row vector

%%% * Box constraint
p_lower = max(v .* B ./ sum(abs(v),2));
p_lower(p_lower <= 0) = 1e-12; % Replace non-positive with a small positive number
p_upper = norm(B, 1) * ones(1, m);
mu_lower = log(p_lower);
mu_upper = log(p_upper);

%%% * Shared Parameters
L_smooth_initial = exp(max(mu_upper)) + (sum(B) / delta_initial); % Initial Lipschitz constant for Adaptive GD/NAG
sigma_calc = min(exp(mu_lower));
if sigma_calc <=0
    warning('Calculated sigma is non-positive. Setting sigma to a small positive value.');
    sigma = 1e-10;
else
    sigma = sigma_calc; % Strong convexity parameter (used by Adaptive NAG)
end
step_size_sub_initial = 1e-3; % Initial step size for Adaptive Subgradient
eta_da = 0.50;               % Step size for Dual Averaging (needs tuning)
gd_step_rule = @(L) 1/L;     % Step size rule for Adaptive GD

%%% * - Initialization
if exist('linear_init_gd', 'file') ~= 2 || exist('linear_init_md', 'file') ~= 2
    error('Initialization functions (linear_init_gd.m or linear_init_md.m) not found.');
end
p0 = linear_init_gd(p_lower,p_upper,sum(B)); % Initial feasible p for dual/subgradient methods
mu0 = log(p0); % Initial mu for adaptive dual methods (GD, NAG)
x0 = linear_init_md(p0,B); % Initial x for primal methods (MD, DA)

disp('--- Running Selected Adaptive/Multi-Turn Algorithms ---');

results = struct(); % Structure to store results

%%% * - solve the problem by Adaptive Vanilla GD
if run_adaptive_gd
    if exist('linear_dual_adaptive_vanilla', 'file') == 2 && exist('linear_dual_gd', 'file') == 2
        disp('Running Adaptive Vanilla GD...');
        adaptive_flag_internal = true; % Enable adaptive checks within linear_dual_gd
        [solution_adpt_gd, time_adpt_gd, iter_adpt_gd, obj_values_adpt_gd, dis_adpt_gd] = ...
            linear_dual_adaptive_vanilla(v, B, mu0, max_iter_per_phase, L_smooth_initial, epsilon, mu_lower, mu_upper, delta_initial, ...
                                         plot_individual_phases, plot_overall_adaptive, p_opt_solver, fval_solver, phase_num, gd_step_rule, adaptive_flag_internal);
        disp(['Adaptive Vanilla GD time: ', num2str(time_adpt_gd), ' seconds, Total Iterations: ', num2str(iter_adpt_gd)]);
        results.adaptive_gd.solution = solution_adpt_gd;
        results.adaptive_gd.obj_values = obj_values_adpt_gd;
        results.adaptive_gd.time = time_adpt_gd;
        results.adaptive_gd.iter = iter_adpt_gd; % Total iterations from all phases
        results.adaptive_gd.distance_final = norm(exp(solution_adpt_gd) - p_opt_solver); % Dual method
    else
        warning('linear_dual_adaptive_vanilla.m or linear_dual_gd.m not found. Skipping Adaptive GD.');
        run_adaptive_gd = false;
    end
end

%%% * - solve the problem by Adaptive NAG
if run_adaptive_nag
    if exist('linear_dual_adaptive', 'file') == 2 && exist('linear_dual_agd', 'file') == 2
        disp('Running Adaptive NAG...');
        adaptive_flag_internal = true; % Enable adaptive checks within linear_dual_agd
        plot_smooth_internal = false; % Disable separate smooth plot from within adaptive call
        [solution_adpt_nag, time_adpt_nag, iter_adpt_nag, obj_values_adpt_nag, dis_adpt_nag] = ...
            linear_dual_adaptive(v, B, mu0, max_iter_per_phase, L_smooth_initial, sigma, epsilon, mu_lower, mu_upper, delta_initial, ...
                                 plot_individual_phases, plot_overall_adaptive, plot_smooth_internal, p_opt_solver, fval_solver, adaptive_flag_internal, phase_num);
        disp(['Adaptive NAG time: ', num2str(time_adpt_nag), ' seconds, Total Iterations: ', num2str(iter_adpt_nag)]);
        results.adaptive_nag.solution = solution_adpt_nag;
        results.adaptive_nag.obj_values = obj_values_adpt_nag;
        results.adaptive_nag.time = time_adpt_nag;
        results.adaptive_nag.iter = iter_adpt_nag; % Total iterations from all phases
        results.adaptive_nag.distance_final = norm(exp(solution_adpt_nag) - p_opt_solver); % Dual method
    else
        warning('linear_dual_adaptive.m or linear_dual_agd.m not found. Skipping Adaptive NAG.');
        run_adaptive_nag = false;
    end
end

%%% * - solve the problem by Adaptive Subgradient
if run_adaptive_sub
    if exist('linear_dual_subgradient_adaptive', 'file') == 2
        disp('Running Adaptive Subgradient...');
        % Note: This function works with p, not mu. Uses p0.
        [solution_adpt_sub, obj_values_adpt_sub, dis_adpt_sub, time_adpt_sub, iter_adpt_sub] = ...
            linear_dual_subgradient_adaptive(v, B, p0, max_iter_adaptive_sub, step_size_sub_initial, epsilon, ...
                                             plot_overall_adaptive, p_opt_solver, fval_solver, switch_step_sub);
        disp(['Adaptive Subgradient time: ', num2str(time_adpt_sub), ' seconds, Iterations: ', num2str(iter_adpt_sub)]);
        results.adaptive_sub.solution = solution_adpt_sub;
        results.adaptive_sub.obj_values = obj_values_adpt_sub;
        results.adaptive_sub.time = time_adpt_sub;
        results.adaptive_sub.iter = iter_adpt_sub; % Total iterations
        results.adaptive_sub.distance_final = norm(solution_adpt_sub - p_opt_solver); % Subgradient returns p
    else
        warning('linear_dual_subgradient_adaptive.m not found. Skipping Adaptive Subgradient.');
        run_adaptive_sub = false;
    end
end

%%% * - solve the problem by Adaptive Mirror Descent
if run_adaptive_md
     if exist('linear_primal_md_adaptive', 'file') == 2
        disp('Running Adaptive Mirror Descent...');
        % Note: This function works with x, initializes with x0.
        [solution_adpt_md, time_adpt_md, iter_adpt_md, obj_values_adpt_md, distance_adpt_md] = ...
            linear_primal_md_adaptive(v, B, x0, epsilon, max_iter_adaptive_md, ...
                                      plot_overall_adaptive, p_opt_solver, fval_solver, switch_step_md);
        disp(['Adaptive MD time: ', num2str(time_adpt_md), ' seconds, Iterations: ', num2str(iter_adpt_md)]);
        results.adaptive_md.solution = solution_adpt_md; % Returns p
        results.adaptive_md.obj_values = obj_values_adpt_md;
        results.adaptive_md.time = time_adpt_md;
        results.adaptive_md.iter = iter_adpt_md; % Total iterations
        results.adaptive_md.distance_final = norm(solution_adpt_md - p_opt_solver); % MD returns p
     else
        warning('linear_primal_md_adaptive.m not found. Skipping Adaptive Mirror Descent.');
        run_adaptive_md = false;
     end
end

%%% * - solve the problem by Dual Averaging (Non-Adaptive)
if run_dual_averaging
     if exist('linear_primal_da', 'file') == 2
        disp('Running Dual Averaging...');
        % Note: Uses x0, non-adaptive iterations
        [solution_da, time_da, iter_da, obj_values_da, distance_da] = ...
            linear_primal_da(v, B, x0, eta_da, epsilon, max_iter_da, ...
                             plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['DA time: ', num2str(time_da), ' seconds, Iterations: ', num2str(iter_da)]);
        results.da.solution = solution_da; % Returns p
        results.da.obj_values = obj_values_da;
        results.da.time = time_da;
        results.da.iter = iter_da; % Total iterations
        results.da.distance_final = norm(solution_da - p_opt_solver); % DA returns p
     else
        warning('linear_primal_da.m not found. Skipping Dual Averaging.');
        run_dual_averaging = false;
     end
end

disp('--- Algorithm Execution Finished ---');

%%% * - plot the combined descent graph (Using Original Objective Gap for comparison)
if run_adaptive_gd || run_adaptive_nag || run_adaptive_sub || run_adaptive_md || run_dual_averaging
    disp('Plotting Objective Gap results...');
    figure; % Create figure for comparison plot
    hold on;
    legend_entries = {}; % Cell array to store legend entries dynamically
    all_obj_values_to_plot = []; % Collect objective values for ylim calculation

    % Plot Adaptive GD results
    if run_adaptive_gd && isfield(results, 'adaptive_gd') && ~isempty(results.adaptive_gd.obj_values)
        total_iters = results.adaptive_gd.iter; % Get total iterations for this method
        semilogy(1:total_iters, abs(results.adaptive_gd.obj_values), '-^', 'DisplayName', 'Adaptive GD', 'LineWidth', 2, 'MarkerSize', 4); % LineWidth=2
        legend_entries{end+1} = 'Adaptive GD';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_gd.obj_values(:))];
    end

    % Plot Adaptive NAG results
    if run_adaptive_nag && isfield(results, 'adaptive_nag') && ~isempty(results.adaptive_nag.obj_values)
        total_iters = results.adaptive_nag.iter;
        semilogy(1:total_iters, abs(results.adaptive_nag.obj_values), '-d', 'DisplayName', 'Adaptive NAG', 'LineWidth', 2, 'MarkerSize', 4); % LineWidth=2
        legend_entries{end+1} = 'Adaptive NAG';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_nag.obj_values(:))];
    end

    % Plot Adaptive Subgradient results
    if run_adaptive_sub && isfield(results, 'adaptive_sub') && ~isempty(results.adaptive_sub.obj_values)
        total_iters = results.adaptive_sub.iter;
        semilogy(1:total_iters, abs(results.adaptive_sub.obj_values), '-s', 'DisplayName', 'Adaptive Subgradient', 'LineWidth', 2, 'MarkerSize', 4); % LineWidth=2
        legend_entries{end+1} = 'Adaptive Subgradient';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_sub.obj_values(:))];
    end

    % Plot Adaptive MD results
    if run_adaptive_md && isfield(results, 'adaptive_md') && ~isempty(results.adaptive_md.obj_values)
        total_iters = results.adaptive_md.iter;
        semilogy(1:total_iters, abs(results.adaptive_md.obj_values), '-o', 'DisplayName', 'Adaptive MD', 'LineWidth', 2, 'MarkerSize', 4); % LineWidth=2
        legend_entries{end+1} = 'Adaptive MD';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_md.obj_values(:))];
    end

    % Plot Dual Averaging results
    if run_dual_averaging && isfield(results, 'da') && ~isempty(results.da.obj_values)
        total_iters = results.da.iter;
        semilogy(1:total_iters, abs(results.da.obj_values), '-*', 'DisplayName', 'Dual Averaging', 'LineWidth', 2, 'MarkerSize', 4); % LineWidth=2
        legend_entries{end+1} = 'Dual Averaging';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.da.obj_values(:))];
    end

    hold off;

    % Customize plot (Matching main_EC_Linear_Synthetic.m style)
    set(gca, 'FontSize', 15); % Axis font size
    xlabel('Iteration', 'FontSize', 20); % X-axis label font size
    ylabel('Original Objective Value Gap F(\cdot) - F*', 'FontSize', 20); % Y-axis label font size
    title_str_obj = sprintf('Multi-Turn Objective Gap Comparison (n=%d, m=%d)', n, m); % Removed delta
    title(title_str_obj, 'FontSize', 20); % Title font size
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'best');
    else
        warning('No algorithms produced results for objective gap plotting.');
    end
    grid on;
    set(gca, 'YScale', 'log');

    % Adjust y-axis limits dynamically based on plotted data
    valid_obj_values = all_obj_values_to_plot(all_obj_values_to_plot > 0 & isfinite(all_obj_values_to_plot));
    if ~isempty(valid_obj_values)
        min_val = max(min(valid_obj_values)*0.1, epsilon/10);
        max_val = max(valid_obj_values)*10;
        if min_val <= 0 || ~isfinite(min_val)
            min_val = epsilon/100; % Fallback small positive value
        end
         if max_val <= min_val || ~isfinite(max_val)
            max_val = max(1, min_val * 100); % Fallback max value
        end
        ylim([min_val, max_val]);
    else
        ylim([epsilon/10, 1]); % Default if no valid data
    end
else
    disp('No algorithms selected or ran successfully for objective gap plot.');
end

% --- No separate smoothed objective plot for multi-turn comparison ---
% The individual adaptive functions might plot internally if flags are set.

%%% * - Calculate and Print final distances to optimal primal solution p*
disp('--- Final Distance to Optimal Primal Solution p* ---');
found_results = false; % Flag to check if any results were printed

if run_adaptive_gd && isfield(results, 'adaptive_gd')
    disp(['Adaptive GD: ', num2str(results.adaptive_gd.distance_final)]);
    found_results = true;
end

if run_adaptive_nag && isfield(results, 'adaptive_nag')
    disp(['Adaptive NAG: ', num2str(results.adaptive_nag.distance_final)]);
    found_results = true;
end

if run_adaptive_sub && isfield(results, 'adaptive_sub')
    disp(['Adaptive Subgradient: ', num2str(results.adaptive_sub.distance_final)]);
    found_results = true;
end

if run_adaptive_md && isfield(results, 'adaptive_md')
    disp(['Adaptive MD: ', num2str(results.adaptive_md.distance_final)]);
    found_results = true;
end

if run_dual_averaging && isfield(results, 'da')
    disp(['Dual Averaging: ', num2str(results.da.distance_final)]);
    found_results = true;
end

if ~found_results % Check if any results were printed
    disp('No results to display.');
end

disp('--- Script Finished ---');
