%%% EC 2025 Submission
%%% ! Version: EC - Linear - Sythetic Data - Multi-Turn/Adaptive Iterative methods
%%% ! Selectable Adaptive Algorithms Comparison
clc
clear

% FOR DUAL AVERAGING, THE OPTIMAL PARAMETERS ARE:
% We only consider the exp generation
% 50*50 0.01
% 100*100 0.005
% 300*300 0.002

% The only thing we need to do is the multiturn experiment

% --- Algorithm Selection ---
run_adaptive_nag = false;  % Set to true to run APM
run_adaptive_sub = false; % Set to true to run Adaptive Tatonnement (was Adaptive Subgradient)
run_adaptive_md = false;  % Set to true to run Adaptive Mirror Descent
run_primal_dual_pdhg = false; % Set to true to run PDHG (was Primal Dual PDHG)
run_primal_md = true; % New: Set to true to run Mirror Descent (was Primal MD)
run_subgradient = false; % New: Set to true to run Tatonnement (was Dual Subgradient)
run_dual_averaging = false; % Set to true to run Dual Averaging (non-adaptive version)
run_adaptive_gd = false;   % Set to true to run Adaptive Vanilla GD

% --------------------------

% --- Load the new dataset ---
v = readmatrix('Dataset/Ratings_kroer.csv') + 0.1;
v = floor(v) + 1; % Preprocess the data
[n, m] = size(v);
B = ones(n, 1); % Unit budget
disp(['Loaded new dataset: Ratings_kroer.csv with n=', num2str(n), ', m=', num2str(m)]);
% --- End new dataset loading ---


% --- Iteration and Tolerance Parameters ---
% General
epsilon = 1e-4; % Stopping criteria based on original objective gap

% Adaptive GD / NAG Parameters
% Todo: check this max_iter_per_phase
max_iter_per_phase = 1000; % Max iterations *per phase* for Adaptive GD/NAG
phase_num = 20;          % Number of phases for Adaptive GD/NAG
delta_initial = 0.1;     % Initial smoothing parameter (used by Adaptive GD/NAG)

% Adaptive Tatonnement / Adaptive Mirror Descent Parameters (was Adaptive Subgradient / MD Parameters)
max_iter_adaptive_sub = 20000; % Max total iterations for Adaptive Tatonnement
max_iter_adaptive_md = 20000;  % Max total iterations for Adaptive Mirror Descent
switch_step_sub = 5000;       % Iteration count to switch step size in Adaptive Tatonnement
switch_step_md = 5000;        % Iteration count to switch step size in Adaptive Mirror Descent

% Dual Averaging Parameters
max_iter_da = 20000;          % Max iterations for Dual Averaging

% New: Mirror Descent Parameters (was Primal MD Parameters)
max_iter_md = 20000;         % Max iterations for Mirror Descent
% Todo: check this eta_md and its adaptive variant
% ! It is correct, but there is please remember to keep the adaptive or not version the same parameters
eta_md = 0.1; % Corrected: Step size parameter for Mirror Descent

% New: Tatonnement Parameters (was Dual Subgradient Parameters)
max_iter_sub = 20000;        % Max iterations for Tatonnement
step_size_sub = 1e-3; % Corrected: Step size parameter for Tatonnement


% PDHG Specific Parameters (Copied from linear_main_wine_syn_singleturn.m)
% Todo: check this for the boosting performance
% Todo: do not let it be too large, otherwise it will be too slow
% ! Still require tuning
max_outer_iter_pdhg = 100; % Number of outer loop restarts for PDHG
K_inner_iter_pdhg = 200; % Inner iterations per outer loop for PDHG
sigma_pdhg = 10/(2*sqrt(n)); % Example initial value, requires tuning
tau_pdhg = 2/(2*sqrt(n)); % Example initial value, requires tuning


% Plotting Flags (Control plotting *within* the called functions, if implemented)
plot_individual_phases = false; % e.g., plot for each GD/NAG phase call
% Todo: check the setting of this parameter
plot_overall_adaptive = false; % e.g., plot overall convergence from within adaptive functions

% --- End Parameters ---

% --- Load pre-computed optimal solver results ---
load('/Users/chjiang/Dropbox/EG_EXP/Linear/WINE_Submission/Efficiency/Solver_Rating.mat', 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');
disp(['Loaded solver results from /Users/chjiang/Dropbox/EG_EXP/Linear/EC_Submission/Efficiency/Solver_Rating.mat']);
disp(['Solver Optimal Value (Original Func): ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
% Todo: Do we need to transpose p_opt_solver?
% p_opt_solver = p_opt_solver'; % Ensure it's a row vector
% --- End pre-computed optimal solver results loading ---

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
step_size_sub_initial = 1e-3; % Initial step size for Adaptive Tatonnement (was Adaptive Subgradient)
eta_da = 0.01;               % Step size for Dual Averaging (needs tuning)
gd_step_rule = @(L) 1/L;     % Step size rule for Adaptive GD

%%% * - Initialization
if exist('linear_init_gd', 'file') ~= 2 || exist('linear_init_md', 'file') ~= 2
    error('Initialization functions (linear_init_gd.m or linear_init_md.m) not found.');
end
p0 = linear_init_gd(p_lower,p_upper,sum(B)); % Initial feasible p for dual/subgradient methods
mu0 = log(p0); % Initial mu for adaptive dual methods (GD, NAG)
x0 = linear_init_md(p0,B); % Initial x for primal methods (MD, DA)

% --- PDHG Specific Initializations ---
x0_for_primal_algs = x0; % Re-using x0 from linear_init_md
p0_for_dual_space_algs = p0; % Re-using p0 from linear_init_gd

% ! Check whether it is standard algorithms for pdhg
t0_for_pdhg = ones(n, 1);
t0_for_pdhg(t0_for_pdhg <= 0) = 1e-12; % Ensure positivity t > 0
y0_for_pdhg = zeros(n, 1); % y can be initialized to zeros
% --- End PDHG Specific Initializations ---


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

%%% * - solve the problem by APM
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

%%% * - solve the problem by Adaptive Tatonnement (was Adaptive Subgradient)
if run_adaptive_sub
    if exist('linear_dual_subgradient_adaptive', 'file') == 2
        disp('Running Adaptive Tatonnement...');
        % Note: This function works with p, not mu. Uses p0.
        [solution_adpt_sub, obj_values_adpt_sub, dis_adpt_sub, time_adpt_sub, iter_adpt_sub] = ...
            linear_dual_subgradient_adaptive(v, B, p0, max_iter_adaptive_sub, step_size_sub_initial, epsilon, ...
                                             plot_overall_adaptive, p_opt_solver, fval_solver, switch_step_sub);
        disp(['Adaptive Tatonnement time: ', num2str(time_adpt_sub), ' seconds, Iterations: ', num2str(iter_adpt_sub)]);
        results.adaptive_sub.solution = solution_adpt_sub;
        results.adaptive_sub.obj_values = obj_values_adpt_sub;
        results.adaptive_sub.time = time_adpt_sub;
        results.adaptive_sub.iter = iter_adpt_sub; % Total iterations
        results.adaptive_sub.distance_final = norm(solution_adpt_sub - p_opt_solver); % Subgradient returns p
    else
        warning('linear_dual_subgradient_adaptive.m not found. Skipping Adaptive Tatonnement.');
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

%%% * - solve the problem by Mirror Descent (was Primal Mirror Descent)
if run_primal_md
    if exist('linear_primal_md', 'file') == 2
        disp('Running Mirror Descent...');
        % Note: Uses x0, non-adaptive iterations. Assumes linear_primal_md returns p.
        plot_flag_md = plot_overall_adaptive; % Use the general plot flag
        [solution_md, time_md, iter_md, obj_values_md, distance_md] = ...
            linear_primal_md(v, B, x0, eta_md, epsilon, max_iter_md, ...
                             plot_flag_md, p_opt_solver, fval_solver);
        disp(['Mirror Descent time: ', num2str(time_md), ' seconds, Iterations: ', num2str(iter_md)]);
        results.primal_md.solution = solution_md; % Returns p
        results.primal_md.obj_values = obj_values_md;
        results.primal_md.time = time_md;
        results.primal_md.iter = iter_md; % Total iterations
        results.primal_md.distance_final = norm(solution_md - p_opt_solver); % MD returns p
    else
        warning('linear_primal_md.m not found. Skipping Mirror Descent.');
        run_primal_md = false;
    end
end

%%% * - solve the problem by Tatonnement (was Dual Subgradient)
if run_subgradient
    if exist('linear_dual_subgradient', 'file') == 2
        disp('Running Tatonnement...');
        % Note: This function works with p, initializes with p0.
        plot_flag_sub = plot_overall_adaptive; % Use the general plot flag
        [solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = ...
            linear_dual_subgradient(v, B, p0, max_iter_sub, step_size_sub, epsilon, ...
                                    plot_flag_sub, p_opt_solver, fval_solver);
        disp(['Tatonnement time: ', num2str(time_sub), ' seconds, Iterations: ', num2str(iter_sub)]);
        results.subgradient.solution = solution_sub; % Subgradient returns p
        results.subgradient.obj_values = obj_values_sub;
        results.subgradient.time = time_sub;
        results.subgradient.iter = iter_sub; % Total iterations
        results.subgradient.distance_final = norm(solution_sub - p_opt_solver); % Compare prices
    else
        warning('linear_dual_subgradient.m not found. Skipping Tatonnement.');
        run_subgradient = false;
    end
end

%%% * - solve the problem by PDHG (was Primal-Dual PDHG)
if run_primal_dual_pdhg
    if exist('linear_primal_dual_pdhg', 'file') == 2
        disp('Running PDHG...');
        % PDHG parameters (x0, t0, p0, y0, max_outer_iter, K_inner_iter, sigma_step, tau_step, epsilon, plot_flag, p_opt_solver, fval_solver)
        [solution_x_pdhg, solution_t_pdhg, solution_p_pdhg, solution_y_pdhg, obj_values_pdhg, dis_p_pdhg, time_pdhg, iter_pdhg] = ...
            linear_primal_dual_pdhg(v, B, x0_for_primal_algs, t0_for_pdhg, p0_for_dual_space_algs, y0_for_pdhg, ...
                                    max_outer_iter_pdhg, K_inner_iter_pdhg, sigma_pdhg, tau_pdhg, epsilon, plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['PDHG time: ', num2str(time_pdhg), ' seconds, Total Inner Iterations: ', num2str(iter_pdhg)]);
        results.pdhg.solution_x = solution_x_pdhg;
        results.pdhg.solution_t = solution_t_pdhg;
        results.pdhg.solution_p = solution_p_pdhg; % Directly from algorithm
        results.pdhg.solution_y = solution_y_pdhg;
        results.pdhg.obj_values = obj_values_pdhg;
        results.pdhg.time = time_pdhg;
        results.pdhg.iter = iter_pdhg;
        results.pdhg.distance_final = norm(solution_p_pdhg - p_opt_solver); % Compare prices
    else
        warning('linear_primal_dual_pdhg.m not found. Skipping PDHG.');
        run_primal_dual_pdhg = false;
    end
end


disp('--- Algorithm Execution Finished ---');

%%% * - plot the combined descent graph (Using Original Objective Gap for comparison)
if run_adaptive_gd || run_adaptive_nag || run_adaptive_sub || run_adaptive_md || run_dual_averaging || run_primal_dual_pdhg || run_primal_md || run_subgradient
    disp('Plotting Objective Gap results...');
    figure; % Create figure for comparison plot
    hold on;
    legend_entries = {}; % Cell array to store legend entries dynamically
    all_obj_values_to_plot = []; % Collect objective values for ylim calculation

    % Get the maximum iteration count across all run algorithms for consistent x-axis
    max_overall_iters = 0;
    if run_adaptive_gd && isfield(results, 'adaptive_gd'), max_overall_iters = max(max_overall_iters, results.adaptive_gd.iter); end
    if run_adaptive_nag && isfield(results, 'adaptive_nag'), max_overall_iters = max(max_overall_iters, results.adaptive_nag.iter); end
    if run_adaptive_sub && isfield(results, 'adaptive_sub'), max_overall_iters = max(max_overall_iters, results.adaptive_sub.iter); end
    if run_adaptive_md && isfield(results, 'adaptive_md'), max_overall_iters = max(max_overall_iters, results.adaptive_md.iter); end
    if run_dual_averaging && isfield(results, 'da'), max_overall_iters = max(max_overall_iters, results.da.iter); end
    if run_primal_md && isfield(results, 'primal_md'), max_overall_iters = max(max_overall_iters, results.primal_md.iter); end
    if run_subgradient && isfield(results, 'subgradient'), max_overall_iters = max(max_overall_iters, results.subgradient.iter); end
    if run_primal_dual_pdhg && isfield(results, 'pdhg'), max_overall_iters = max(max_overall_iters, results.pdhg.iter); end


    % Plot Adaptive GD results
    if run_adaptive_gd && isfield(results, 'adaptive_gd') && ~isempty(results.adaptive_gd.obj_values)
        total_iters = results.adaptive_gd.iter; % Get total iterations for this method
        plot_range_algo = 1:(total_iters - 1); % Apply -1 here
        % Ensure plot_range_algo is not empty or invalid if total_iters is 1
        if ~isempty(plot_range_algo) && total_iters > 0
            h = semilogy(plot_range_algo, abs(results.adaptive_gd.obj_values(plot_range_algo)), '-^', 'DisplayName', 'Adaptive GD', 'LineWidth', 2, 'MarkerSize', 4);
            legend_entries{end+1} = 'Adaptive GD';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_gd.obj_values(:))];
        end
    end

    % Plot Adaptive NAG results
    if run_adaptive_nag && isfield(results, 'adaptive_nag') && ~isempty(results.adaptive_nag.obj_values)
        total_iters = results.adaptive_nag.iter;
        plot_range_algo = 1:(total_iters - 1); % Apply -1 here
        if ~isempty(plot_range_algo) && total_iters > 0
            h = semilogy(plot_range_algo, abs(results.adaptive_nag.obj_values(plot_range_algo)), '-d', 'DisplayName', 'Adaptive NAG', 'LineWidth', 2, 'MarkerSize', 4);
            legend_entries{end+1} = 'APM';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_nag.obj_values(:))];
        end
    end

    % Plot Adaptive Tatonnement results (was Adaptive Subgradient)
    if run_adaptive_sub && isfield(results, 'adaptive_sub') && ~isempty(results.adaptive_sub.obj_values)
        total_iters = results.adaptive_sub.iter;
        plot_range_algo = 1:(total_iters - 1); % Apply -1 here
        if ~isempty(plot_range_algo) && total_iters > 0
            h = semilogy(plot_range_algo, abs(results.adaptive_sub.obj_values(plot_range_algo)), '-s', 'DisplayName', 'Adaptive Tatonnement', 'LineWidth', 2, 'MarkerSize', 4);
            legend_entries{end+1} = 'Adaptive Tatonnement';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_sub.obj_values(:))];
        end
    end

    % Plot Adaptive MD results
    if run_adaptive_md && isfield(results, 'adaptive_md') && ~isempty(results.adaptive_md.obj_values)
        total_iters = results.adaptive_md.iter;
        plot_range_algo = 1:(total_iters - 1); % Apply -1 here
        if ~isempty(plot_range_algo) && total_iters > 0
            h = semilogy(plot_range_algo, abs(results.adaptive_md.obj_values(plot_range_algo)), '-o', 'DisplayName', 'Adaptive Mirror Descent', 'LineWidth', 2, 'MarkerSize', 4);
            legend_entries{end+1} = 'Adaptive MD';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_md.obj_values(:))];
        end
    end

    % Plot Dual Averaging results
    if run_dual_averaging && isfield(results, 'da') && ~isempty(results.da.obj_values)
        total_iters = results.da.iter;
        plot_range_algo = 1:(total_iters - 1); % Apply -1 here
        if ~isempty(plot_range_algo) && total_iters > 0
            h = semilogy(plot_range_algo, abs(results.da.obj_values(plot_range_algo)), '-*', 'DisplayName', 'Dual Averaging', 'LineWidth', 2, 'MarkerSize', 4);
            legend_entries{end+1} = 'Dual Averaging';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.da.obj_values(:))];
        end
    end

    % Plot Mirror Descent results (was Primal MD)
    if run_primal_md && isfield(results, 'primal_md') && ~isempty(results.primal_md.obj_values)
        total_iters = results.primal_md.iter;
        plot_range_algo = 1:(total_iters - 1);
        if ~isempty(plot_range_algo) && total_iters > 0
            h = semilogy(plot_range_algo, abs(results.primal_md.obj_values(plot_range_algo)), '-x', 'DisplayName', 'Mirror Descent', 'LineWidth', 2, 'MarkerSize', 4);
            legend_entries{end+1} = 'Mirror Descent';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.primal_md.obj_values(:))];
        end
    end

    % Plot Tatonnement results (was Dual Subgradient)
    if run_subgradient && isfield(results, 'subgradient') && ~isempty(results.subgradient.obj_values)
        total_iters = results.subgradient.iter;
        plot_range_algo = 1:(total_iters - 1);
        if ~isempty(plot_range_algo) && total_iters > 0
            h = semilogy(plot_range_algo, abs(results.subgradient.obj_values(plot_range_algo)), '-v', 'DisplayName', 'Tatonnement', 'LineWidth', 2, 'MarkerSize', 4);
            legend_entries{end+1} = 'Tatonnement';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.subgradient.obj_values(:))];
        end
    end

    % Plot PDHG results (was Primal-Dual PDHG)
    if run_primal_dual_pdhg && isfield(results, 'pdhg') && ~isempty(results.pdhg.obj_values)
        total_iters = results.pdhg.iter;
        plot_range_algo = 1:(total_iters - 1); % Apply -1 here
        if ~isempty(plot_range_algo) && total_iters > 0
            h = semilogy(plot_range_algo, abs(results.pdhg.obj_values(plot_range_algo)), '-<', 'DisplayName', 'PDHG', 'LineWidth', 2, 'MarkerSize', 4);
            legend_entries{end+1} = 'PDHG';
            all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.pdhg.obj_values(:))];
        end
    end

    hold off;

    % Customize plot (Matching singleturn style)
    set(gca, 'FontSize', 15); % Axis font size
    xlabel('Iteration', 'FontSize', 20); % X-axis label font size
    ylabel('Objective Value Gap', 'FontSize', 20); % Y-axis label font size
    title_str_obj = sprintf('n=%d, m=%d', n, m);
    title(title_str_obj, 'FontSize', 25); % Title font size
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'northeast', 'FontSize', 10); %Todo: used to be best+15
    else
        warning('No algorithms produced results for objective gap plotting.');
    end
    grid on;
    set(gca, 'YScale', 'log','LineWidth', 2); %Todo: check the revision here
    % Determine dynamic x-axis limit
    max_plot_iters = 0;
    if run_adaptive_gd && isfield(results, 'adaptive_gd'), max_plot_iters = max(max_plot_iters, results.adaptive_gd.iter - 1); end
    if run_adaptive_nag && isfield(results, 'adaptive_nag'), max_plot_iters = max(max_plot_iters, results.adaptive_nag.iter - 1); end
    if run_adaptive_sub && isfield(results, 'adaptive_sub'), max_plot_iters = max(max_plot_iters, results.adaptive_sub.iter - 1); end
    if run_adaptive_md && isfield(results, 'adaptive_md'), max_plot_iters = max(max_plot_iters, results.adaptive_md.iter - 1); end
    if run_dual_averaging && isfield(results, 'da'), max_plot_iters = max(max_plot_iters, results.da.iter - 1); end
    if run_primal_md && isfield(results, 'primal_md'), max_plot_iters = max(max_plot_iters, results.primal_md.iter - 1); end
    if run_subgradient && isfield(results, 'subgradient'), max_plot_iters = max(max_plot_iters, results.subgradient.iter - 1); end
    if run_primal_dual_pdhg && isfield(results, 'pdhg'), max_plot_iters = max(max_plot_iters, results.pdhg.iter - 1); end
    
    if max_plot_iters > 0
        xlim([1, max_plot_iters]); % Set x-axis limit to max iterations
    else
        xlim([1, 100]); % Default if no algorithms ran or less than 1 iteration
    end

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
    disp(['APM: ', num2str(results.adaptive_nag.distance_final)]);
    found_results = true;
end

if run_adaptive_sub && isfield(results, 'adaptive_sub')
    disp(['Adaptive Tatonnement: ', num2str(results.adaptive_sub.distance_final)]);
    found_results = true;
end

if run_adaptive_md && isfield(results, 'adaptive_md')
    disp(['Adaptive Mirror Descent: ', num2str(results.adaptive_md.distance_final)]);
    found_results = true;
end

if run_dual_averaging && isfield(results, 'da')
    disp(['Dual Averaging: ', num2str(results.da.distance_final)]);
    found_results = true;
end

if run_primal_md && isfield(results, 'primal_md')
    disp(['Mirror Descent: ', num2str(results.primal_md.distance_final)]);
    found_results = true;
end

if run_subgradient && isfield(results, 'subgradient')
    disp(['Tatonnement: ', num2str(results.subgradient.distance_final)]);
    found_results = true;
end

if run_primal_dual_pdhg && isfield(results, 'pdhg')
    disp(['PDHG: ', num2str(results.pdhg.distance_final)]);
    found_results = true;
end

if ~found_results % Check if any results were printed
    disp('No results to display.');
end

disp('--- Script Finished ---');