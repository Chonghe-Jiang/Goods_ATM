%%% EC 2025 Submission
%%% ! Version: EC - Linear - Sythetic Data - Iterative method
%%% ! Optimal Version for Submission
%%% ! Includes Adaptive Vanilla Gradient Descent (GD) comparison
clc
clear

% Define problem parameters
n = 50;  % Number of rows
m = 50;   % Number of columns
B = ones(n,1);

% Define the folder name
dataset_folder = 'synthetica_dataset';

% Iteration and tolerance parameters
max_iter_sub = 20000; % Max iterations for Subgradient
max_iter_md = 20000;  % Max iterations for Mirror Descent
epsilon = 1e-4; % Stopping criteria based on original objective gap
plot_flag = false; % Set to true to see individual algorithm plots (like each GD phase)
adaptive_plot_flag = true; % Set to true to see the overall plot for the adaptive method

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
    [p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
    save(solver_filepath, 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');  % Save solver results
    disp(['Solved and saved solver results to ', solver_filepath]);
end
disp(['Solver Optimal Value (Original Func): ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
p_opt_solver = p_opt_solver'; % Ensure it's a row vector

%%% * Box constraint
p_lower = max(v .* B ./ sum(abs(v),2));
p_upper = norm(B, 1) * ones(1, m);
mu_lower = log(p_lower);
mu_upper = log(p_upper);

%%% * Parameters
delta_initial = 0.1; % Initial smoothing parameter for adaptive method
L_smooth_initial = exp(max(mu_upper)) + (sum(B) / delta_initial); % Initial Lipschitz constant for smoothed gradient
step_size_sub = 1e-3; % Subgradient stepsize (needs tuning)
eta = 1;  % Mirror Descent stepsize (needs tuning)
% Define step size rule for adaptive vanilla GD (e.g., 1/L)
gd_step_rule = @(L) 1/L;
% Alternative: fixed step size (might need tuning)
% gd_step_rule = 1e-4;

%%% * - Initialization
p0 = linear_init_gd(p_lower,p_upper,sum(B)); % Initial feasible p
mu0 = log(p0); % Initial mu
x0 = linear_init_md(p0,B); % Initial x for Mirror Descent (assuming linear_init_md.m exists)

%%% * - solve the problem by subgradient (Tatonnement) - Operates on Original Objective
[solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = linear_dual_subgradient(v, B, p0, max_iter_sub, step_size_sub, epsilon, plot_flag, p_opt_solver, fval_solver);
disp(['Subgradient time: ', num2str(time_sub), ' seconds, Iterations: ', num2str(iter_sub)]);

%%% * - solve the problem by mirror descent - Primal Method
[solution_md, time_md, iter_md, obj_values_md, distance_md] = linear_primal_md(v, B, x0, eta, epsilon, max_iter_md, plot_flag, p_opt_solver, fval_solver);
disp(['MD time: ', num2str(time_md), ' seconds, Iterations: ', num2str(iter_md)]);

% %%% * - solve the problem by Vanilla Gradient Descent (Direct Call - Now replaced by Adaptive version)
% max_iter_gd = 20000;
% step_size_gd = 1 / L_smooth_initial;
% [solution_gd, time_gd, iter_gd, obj_values_gd, f_smooth_values_gd, dis_gd, convergence_gd] = linear_dual_gd(v, B, mu0, max_iter_gd, L_smooth_initial, epsilon, mu_lower, mu_upper, delta_initial, plot_flag, p_opt_solver, fval_solver, step_size_gd);
% disp(['Vanilla GD (Direct) time: ', num2str(time_gd), ' seconds, Iterations: ', num2str(iter_gd)]);

% %%% * - solve the problem by Adaptive AGD (ATM) - Commented out to focus on Vanilla Adaptive
max_iter_adaptive_agd = 4500; % Max iterations for adaptive AGD (original)
sigma = min(exp(mu_lower)); % Strong convexity parameter (used by AGD)
adaptive = true; % share the same adaptive parameter
phase_num_agd = 20;
[solution_adaptive_agd, total_time_adaptive_agd, total_iter_adaptive_agd, obj_values_adaptive_agd, dis_adaptive_agd] = linear_dual_adaptive(v, B, mu0, max_iter_adaptive_agd, L_smooth_initial, sigma, epsilon, mu_lower, mu_upper, delta_initial, plot_flag, adaptive_plot_flag, false, p_opt_solver, fval_solver, adaptive, phase_num_agd);
disp(['Adaptive AGD (ATM) time: ', num2str(total_time_adaptive_agd), ' seconds, Iterations: ', num2str(total_iter_adaptive_agd)]);

%%% * - solve the problem by Adaptive Vanilla Gradient Descent (NEW)
max_iter_per_phase_adaptive_vanilla = 4500; % Max iterations *per phase* for Adaptive Vanilla GD
phase_num_adaptive_vanilla = 20; % Number of phases for Adaptive Vanilla GD
[solution_adaptive_vanilla, total_time_adaptive_vanilla, total_iter_adaptive_vanilla, obj_values_adaptive_vanilla, dis_adaptive_vanilla] = ...
    linear_dual_adaptive_vanilla(v, B, mu0, max_iter_per_phase_adaptive_vanilla, L_smooth_initial, epsilon, mu_lower, mu_upper, delta_initial, ...
                                 plot_flag, adaptive_plot_flag, p_opt_solver, fval_solver, phase_num_adaptive_vanilla, gd_step_rule, adaptive);
disp(['Adaptive Vanilla GD time: ', num2str(total_time_adaptive_vanilla), ' seconds, Iterations: ', num2str(total_iter_adaptive_vanilla)]);


%%% * - plot the combined descent graph (Using Original Objective Gap for comparison)
figure;
hold on;

% Plot Mirror Descent (note: obj_values_md might be different, e.g., primal gap)
semilogy(1:length(obj_values_md), abs(obj_values_md), '-o', 'DisplayName', 'Mirror Descent (Primal Gap?)', 'LineWidth', 1.5, 'MarkerSize', 4);

% Plot Subgradient (Tatonnement)
semilogy(1:length(obj_values_sub), abs(obj_values_sub), '-s', 'DisplayName', 'Tatonnement (Orig. Obj. Gap)', 'LineWidth', 1.5, 'MarkerSize', 4);

% Plot Adaptive Vanilla Gradient Descent (NEW)
semilogy(1:length(obj_values_adaptive_vanilla), abs(obj_values_adaptive_vanilla), '-^', 'DisplayName', 'Adaptive Vanilla GD (Orig. Obj. Gap)', 'LineWidth', 1.5, 'MarkerSize', 4);

% % Plot Adaptive AGD (ATM) - If uncommented above
if exist('obj_values_adaptive_agd', 'var')
   semilogy(1:length(obj_values_adaptive_agd), abs(obj_values_adaptive_agd), '-d', 'DisplayName', 'ATM (Orig. Obj. Gap)', 'LineWidth', 1.5, 'MarkerSize', 4);
end

hold off;

% Customize plot
set(gca, 'FontSize', 12);
xlabel('Iteration', 'FontSize', 14);
ylabel('Absolute Original Objective Value Gap |f(\cdot) - f*|', 'FontSize', 14); % Clarified label
title(['Convergence Comparison (n=', num2str(n), ', m=', num2str(m), ')'], 'FontSize', 16);
legend show;
grid on;
set(gca, 'YScale', 'log');
% Adjust y-axis limits dynamically
all_obj_values = [abs(obj_values_md(:)); abs(obj_values_sub(:)); abs(obj_values_adaptive_vanilla(:))];
% if exist('obj_values_adaptive_agd', 'var')
%     all_obj_values = [all_obj_values; abs(obj_values_adaptive_agd(:))];
% end
valid_obj_values = all_obj_values(all_obj_values > 0 & isfinite(all_obj_values)); % Exclude zero/NaN/Inf
if ~isempty(valid_obj_values)
    ylim([max(min(valid_obj_values)*0.1, epsilon/10), max(valid_obj_values)*10]);
else
    ylim([epsilon/10, 1]); % Default if no valid data
end


%%% * - Calculate final distances to optimal primal solution
distance_md_to_solver = norm(solution_md - p_opt_solver); % MD solution is primal
distance_sub_to_solver = norm(solution_sub - p_opt_solver); % Subgradient solution is primal (p)
distance_adaptive_vanilla_to_solver = norm(exp(solution_adaptive_vanilla) - p_opt_solver); % Adaptive Vanilla GD solution is dual (mu)
% if exist('solution_adaptive_agd', 'var')
%     distance_adaptive_agd_to_solver = norm(exp(solution_adaptive_agd) - p_opt_solver); % Adaptive AGD solution is dual (mu)
% end

% Print the distances
disp('--- Final Distance to Optimal Primal Solution p* ---');
disp(['Mirror Descent: ', num2str(distance_md_to_solver)]);
disp(['Tatonnement (Subgradient): ', num2str(distance_sub_to_solver)]);
disp(['Adaptive Vanilla GD: ', num2str(distance_adaptive_vanilla_to_solver)]);
% if exist('distance_adaptive_agd_to_solver', 'var')
%     disp(['ATM (Adaptive AGD): ', num2str(distance_adaptive_agd_to_solver)]);
% end