%%% EC 2025 Submission
%%% ! Version: EC - Linear - Sythetic Data - Iterative method
%%% ! Optimal Version for Submission
%%% ! Includes Direct GD, NAG, MD, DA, Subgradient comparison with selectable algorithms
%%% ! Added integrated smoothed objective plot when plot_flag_smooth = true
clc
clear

% --- Algorithm Selection ---
run_subgradient = false;      % Set to true to run Subgradient (Tatonnement)
run_mirror_descent = false;   % Set to true to run Mirror Descent (Primal)
run_dual_averaging = true;    % Set to true to run Dual Averaging (Primal) % <-- ADDED
run_gd = false;              % Set to true to run Vanilla GD (Dual Smoothed)
run_agd = true;             % Set to true to run NAG (Dual Smoothed)
% --------------------------

% Define problem parameters
n = 50;  % Number of rows
m = 50;   % Number of columns
B = ones(n,1);

% Define the folder name
dataset_folder = 'synthetica_dataset';

% Iteration and tolerance parameters
max_iter_sub = 20000; % Max iterations for Subgradient
max_iter_md = 20000;  % Max iterations for Mirror Descent
max_iter_da = 20000;  % Max iterations for Dual Averaging % <-- ADDED
max_iter_gd = 20000;  % Max iterations for Vanilla GD
max_iter_agd = 20000; % Max iterations for NAG (adjust as needed)
epsilon = 1e-1; % Stopping criteria based on original objective gap
delta = 0.01; % Smoothing parameter (used by GD and NAG)
plot_flag = false; % Set to true to see individual algorithm plots (Original Gap, Smooth Value, Distance)
plot_flag_smooth = true; % Set to true to see separate plot for smoothed objective value comparison (for GD and NAG)
% adaptive_plot_flag is not used in this version

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
    % Solve the problem using the solver (assuming linear_dual_solver.m exists)
    % Make sure linear_dual_solver is available in the path or current folder
    if exist('linear_dual_solver', 'file') == 2
        [p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
        save(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Save solver results
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
% Handle cases where p_lower might be zero or negative before taking log
p_lower(p_lower <= 0) = 1e-12; % Replace non-positive with a small positive number
p_upper = norm(B, 1) * ones(1, m);
mu_lower = log(p_lower);
mu_upper = log(p_upper);

%%% * Parameters
L_smooth = exp(max(mu_upper)) + (sum(B) / delta); % Lipschitz constant for smoothed gradient
step_size_sub = 1e-3; % Subgradient stepsize (needs tuning)
eta_md = 1;  % Mirror Descent stepsize (needs tuning)
eta_da = 0.50;  % Dual Averaging stepsize (needs tuning, might differ from eta_md) % <-- ADDED
step_size_gd = 1 / L_smooth; % Step size for Vanilla GD (can be tuned)
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
% Ensure linear_init_gd and linear_init_md are available
if exist('linear_init_gd', 'file') ~= 2
    error('linear_init_gd.m not found.');
end
if exist('linear_init_md', 'file') ~= 2
    error('linear_init_md.m not found.');
end
p0 = linear_init_gd(p_lower,p_upper,sum(B)); % Initial feasible p for dual methods
mu0 = log(p0); % Initial mu for dual methods (GD, NAG)
x0 = linear_init_md(p0,B); % Initial x for primal methods (MD, DA)

disp('--- Running Selected Algorithms ---');

results = struct(); % Structure to store results of run algorithms

%%% * - solve the problem by subgradient (Tatonnement) - Operates on Original Objective
if run_subgradient
    if exist('linear_dual_subgradient', 'file') == 2
        disp('Running Subgradient (Tatonnement)...');
        [solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = linear_dual_subgradient(v, B, p0, max_iter_sub, step_size_sub, epsilon, plot_flag, p_opt_solver, fval_solver);
        disp(['Subgradient time: ', num2str(time_sub), ' seconds, Iterations: ', num2str(iter_sub)]);
        results.subgradient.solution = solution_sub;
        results.subgradient.obj_values = obj_values_sub;
        results.subgradient.time = time_sub;
        results.subgradient.iter = iter_sub;
        results.subgradient.distance_final = norm(exp(solution_sub) - p_opt_solver); % Subgradient works on dual (mu), so exp() needed
    else
        warning('linear_dual_subgradient.m not found. Skipping Subgradient.');
        run_subgradient = false; % Disable plotting/reporting for this algorithm
    end
end

%%% * - solve the problem by mirror descent - Primal Method
if run_mirror_descent
     if exist('linear_primal_md', 'file') == 2
        disp('Running Mirror Descent...');
        [solution_md, time_md, iter_md, obj_values_md, distance_md] = linear_primal_md(v, B, x0, eta_md, epsilon, max_iter_md, plot_flag, p_opt_solver, fval_solver);
        disp(['MD time: ', num2str(time_md), ' seconds, Iterations: ', num2str(iter_md)]);
        results.md.solution = solution_md;
        results.md.obj_values = obj_values_md;
        results.md.time = time_md;
        results.md.iter = iter_md;
        results.md.distance_final = norm(solution_md - p_opt_solver); % MD solution is primal p
     else
        warning('linear_primal_md.m not found. Skipping Mirror Descent.');
        run_mirror_descent = false;
     end
end

%%% * - solve the problem by dual averaging - Primal Method % <-- ADDED SECTION
if run_dual_averaging
     if exist('linear_primal_da', 'file') == 2
        disp('Running Dual Averaging...');
        % Note: Using eta_da parameter
        [solution_da, time_da, iter_da, obj_values_da, distance_da] = linear_primal_da(v, B, x0, eta_da, epsilon, max_iter_da, plot_flag, p_opt_solver, fval_solver);
        disp(['DA time: ', num2str(time_da), ' seconds, Iterations: ', num2str(iter_da)]);
        results.da.solution = solution_da;
        results.da.obj_values = obj_values_da;
        results.da.time = time_da;
        results.da.iter = iter_da;
        results.da.distance_final = norm(solution_da - p_opt_solver); % DA solution is primal p
     else
        warning('linear_primal_da.m not found. Skipping Dual Averaging.');
        run_dual_averaging = false;
     end
end

%%% * - solve the problem by Vanilla Gradient Descent (Direct Call) - Dual Smoothed
if run_gd
    if exist('linear_dual_gd', 'file') == 2
        disp('Running Vanilla GD...');
        % Call GD, receive f_smooth_values_gd
        [solution_gd, time_gd, iter_gd, obj_values_gd, f_smooth_values_gd, dis_gd, convergence_gd] = ...
            linear_dual_gd(v, B, mu0, max_iter_gd, L_smooth, epsilon, mu_lower, mu_upper, delta, ...
                           plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, step_size_gd, adaptive);
        disp(['Vanilla GD time: ', num2str(time_gd), ' seconds, Iterations: ', num2str(iter_gd)]);
        results.gd.solution = solution_gd;
        results.gd.obj_values = obj_values_gd;
        results.gd.f_smooth_values = f_smooth_values_gd; % Store smoothed values
        results.gd.time = time_gd;
        results.gd.iter = iter_gd;
        results.gd.distance_final = norm(exp(solution_gd) - p_opt_solver); % GD solution is dual mu
    else
        warning('linear_dual_gd.m not found. Skipping Vanilla GD.');
        run_gd = false;
    end
end

%%% * - solve the problem by Nesterov's Accelerated Gradient (NAG - Direct Call) - Dual Smoothed
if run_agd % Keep variable name run_agd for simplicity, but logic refers to NAG now
     if exist('linear_dual_agd', 'file') == 2 % Assuming the file name remains linear_dual_agd.m
        disp('Running NAG...');
        % Call AGD/NAG, receive f_smooth_values_agd (ASSUMING linear_dual_agd was modified)
        [solution_agd, time_agd, iter_agd, obj_values_agd, f_smooth_values_agd, dis_agd, convergence_agd] = ...
            linear_dual_agd(v, B, mu0, max_iter_agd, L_smooth, sigma, epsilon, mu_lower, mu_upper, delta, ...
                            plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive);
        disp(['NAG time: ', num2str(time_agd), ' seconds, Iterations: ', num2str(iter_agd)]); % Changed output text
        results.agd.solution = solution_agd; % Keep internal struct name agd
        results.agd.obj_values = obj_values_agd;
        results.agd.f_smooth_values = f_smooth_values_agd; % Store smoothed values
        results.agd.time = time_agd;
        results.agd.iter = iter_agd;
        results.agd.distance_final = norm(exp(solution_agd) - p_opt_solver); % NAG solution is dual mu
     else
        warning('linear_dual_agd.m not found or not modified to return smoothed values. Skipping NAG.');
        run_agd = false;
     end
end

disp('--- Algorithm Execution Finished ---');

%%% * - plot the combined descent graph (Using Original Objective Gap for comparison)
if run_subgradient || run_mirror_descent || run_dual_averaging || run_gd || run_agd % Only plot if at least one algorithm ran % <-- ADDED run_dual_averaging
    disp('Plotting Objective Gap results...');
    figure; % Create figure for comparison plot
    hold on;
    legend_entries = {}; % Cell array to store legend entries dynamically
    all_obj_values_to_plot = []; % Collect objective values for ylim calculation

    if run_mirror_descent && isfield(results, 'md') && ~isempty(results.md.obj_values)
        semilogy(1:length(results.md.obj_values), abs(results.md.obj_values), '-o', 'DisplayName', 'Mirror Descent (Primal Gap)', 'LineWidth', 1.5, 'MarkerSize', 4); % Updated label
        legend_entries{end+1} = 'Mirror Descent (Primal Gap)';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.md.obj_values(:))];
    end

    if run_dual_averaging && isfield(results, 'da') && ~isempty(results.da.obj_values) % <-- ADDED DA PLOTTING
        semilogy(1:length(results.da.obj_values), abs(results.da.obj_values), '-*', 'DisplayName', 'Dual Averaging (Primal Gap)', 'LineWidth', 1.5, 'MarkerSize', 4);
        legend_entries{end+1} = 'Dual Averaging (Primal Gap)';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.da.obj_values(:))];
    end

    if run_subgradient && isfield(results, 'subgradient') && ~isempty(results.subgradient.obj_values)
        semilogy(1:length(results.subgradient.obj_values), abs(results.subgradient.obj_values), '-s', 'DisplayName', 'Tatonnement (Dual Gap)', 'LineWidth', 1.5, 'MarkerSize', 4); % Updated label
        legend_entries{end+1} = 'Tatonnement (Dual Gap)';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.subgradient.obj_values(:))];
    end

    if run_gd && isfield(results, 'gd') && ~isempty(results.gd.obj_values)
        semilogy(1:length(results.gd.obj_values), abs(results.gd.obj_values), '-^', 'DisplayName', 'Vanilla GD (Dual Gap)', 'LineWidth', 1.5, 'MarkerSize', 4); % Updated label
        legend_entries{end+1} = 'Vanilla GD (Dual Gap)';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.gd.obj_values(:))];
    end

    if run_agd && isfield(results, 'agd') && ~isempty(results.agd.obj_values)
        semilogy(1:length(results.agd.obj_values), abs(results.agd.obj_values), '-d', 'DisplayName', 'NAG (Dual Gap)', 'LineWidth', 1.5, 'MarkerSize', 4); % Updated label
        legend_entries{end+1} = 'NAG (Dual Gap)';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.agd.obj_values(:))];
    end

    hold off;

    % Customize plot
    set(gca, 'FontSize', 12);
    xlabel('Iteration', 'FontSize', 14);
    % Requirement 4: Remove "Absolute" from Y-axis label
    ylabel('Original Objective Value Gap |F(\cdot) - F*|', 'FontSize', 14); % Capital F, added absolute value bars
    % Requirement 5: Remove delta from title
    title_str_obj = sprintf('Objective Gap Comparison (n=%d, m=%d)', n, m);
    title(title_str_obj, 'FontSize', 16);
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'best'); % Show legend only for plotted algorithms
    else
        warning('No algorithms produced results for objective gap plotting.');
    end
    grid on;
    set(gca, 'YScale', 'log');

    % Adjust y-axis limits dynamically based on plotted data
    valid_obj_values = all_obj_values_to_plot(all_obj_values_to_plot > 0 & isfinite(all_obj_values_to_plot)); % Exclude zero/NaN/Inf
    if ~isempty(valid_obj_values)
        min_val = max(min(valid_obj_values)*0.1, epsilon/10); % Ensure min_val is positive
        max_val = max(valid_obj_values)*10;
         % Ensure min_val < max_val before setting ylim
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

%%% * -- Plot the integrated smoothed objective function graph -- %%%
% No changes needed here as DA does not use the smoothed objective
if plot_flag_smooth && run_gd && run_agd % Check if flag is true AND both GD and NAG ran
    % Check if smoothed results are available for both
    if isfield(results, 'gd') && isfield(results.gd, 'f_smooth_values') && ~isempty(results.gd.f_smooth_values) && ...
       isfield(results, 'agd') && isfield(results.agd, 'f_smooth_values') && ~isempty(results.agd.f_smooth_values)

        disp('Plotting Smoothed Objective Function Comparison...');
        figure; % Create a new figure specifically for this plot
        hold on;

        % Plot GD Smoothed Values
        plot(1:length(results.gd.f_smooth_values), results.gd.f_smooth_values, '-^', 'DisplayName', 'Vanilla GD', 'LineWidth', 1.5, 'MarkerSize', 4);

        % Plot NAG Smoothed Values
        plot(1:length(results.agd.f_smooth_values), results.agd.f_smooth_values, '-d', 'DisplayName', 'NAG', 'LineWidth', 1.5, 'MarkerSize', 4); % Changed DisplayName

        hold off;

        % Customize plot
        set(gca, 'FontSize', 12);
        xlabel('Iteration', 'FontSize', 14);
        % Requirement 2: Capital F in Y-axis label
        ylabel('Smoothed Objective Value F_{\delta}(\mu)', 'FontSize', 14);
        % Requirement 1 & 5: New Title without delta
        title_str_smooth = sprintf('Comparison on Problem (P_{\\delta}) (n=%d, m=%d)', n, m);
        title(title_str_smooth, 'FontSize', 16);
        legend('Location', 'best');
        grid on;
        % Optional: Adjust Y-axis limits if needed
        % all_smooth_vals = [results.gd.f_smooth_values(:); results.agd.f_smooth_values(:)];
        % ylim([min(all_smooth_vals)*0.99, max(all_smooth_vals)*1.01]);

    else
        warning('plot_flag_smooth is true, but smoothed objective results for both GD and NAG are not available. Skipping plot.');
    end
elseif plot_flag_smooth
    disp('plot_flag_smooth is true, but both GD and NAG were not selected or did not run successfully. Skipping smoothed objective plot.');
end


%%% * - Calculate and Print final distances to optimal primal solution p*
disp('--- Final Distance to Optimal Primal Solution p* ---');
found_results = false; % Flag to check if any results were printed

if run_mirror_descent && isfield(results, 'md')
    disp(['Mirror Descent: ', num2str(results.md.distance_final)]);
    found_results = true;
end

if run_dual_averaging && isfield(results, 'da') % <-- ADDED DA REPORTING
    disp(['Dual Averaging: ', num2str(results.da.distance_final)]);
    found_results = true;
end

if run_subgradient && isfield(results, 'subgradient')
    disp(['Tatonnement (Subgradient): ', num2str(results.subgradient.distance_final)]);
    found_results = true;
end

if run_gd && isfield(results, 'gd')
    disp(['Vanilla GD: ', num2str(results.gd.distance_final)]);
    found_results = true;
end

if run_agd && isfield(results, 'agd')
    disp(['NAG: ', num2str(results.agd.distance_final)]); % Changed display text
    found_results = true;
end

if ~found_results % Check if any results were printed
    disp('No results to display.');
end

disp('--- Script Finished ---');

