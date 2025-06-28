%%% EC 2025 Submission
%%% ! Version: EC - Linear - Sythetic Data - Multi-Turn/Adaptive Iterative methods
%%% ! Selectable Adaptive Algorithms Comparison
%%% ! MODIFIED to include Adaptive PDHG and full algorithm suite
clc
clear
close all;

% --- Algorithm Selection ---
% Set flags to true to run the corresponding algorithms.
run_adaptive_nag = true;
run_adaptive_sub = true;
run_adaptive_md = true;
run_primal_dual_pdhg_fixed = false; % Set to true to run the original FIXED step-size PDHG
run_primal_dual_pdhg_adaptive = false; % Set to true to run the NEW ADAPTIVE step-size PDHG
run_primal_md = true;
run_subgradient = true;
run_dual_averaging = true;
run_adaptive_gd = false;

% --- Load the new dataset ---
v = readmatrix('Dataset/Ratings_kroer.csv') + 0.1;
v = floor(v) + 1; % Preprocess the data
[n, m] = size(v);
B = ones(n, 1); % Unit budget
disp(['Loaded new dataset: Ratings_kroer.csv with n=', num2str(n), ', m=', num2str(m)]);

% --- Iteration and Tolerance Parameters ---
% General
epsilon = 1e-4; % Stopping criteria based on original objective gap

% PDHG Specific Parameters
max_outer_iter_pdhg = 200; % Number of outer loop restarts for PDHG
K_inner_iter_pdhg = 100;   % Inner iterations per outer loop for PDHG

% --- Parameters for Fixed Step-Size PDHG ---
sigma_pdhg_fixed = 1/(sqrt(n)); % Example fixed value, requires tuning
tau_pdhg_fixed = 1/(sqrt(n));   % Example fixed value, requires tuning

% --- Parameters for Adaptive Step-Size PDHG (as per paper) ---
eta_initial_pdhg_adaptive = 1/(sqrt(n));    % Initial overall step-size
omega_initial_pdhg_adaptive = 1.0;  % Initial primal-dual weight (1.0 is balanced)
theta_pdhg_adaptive = 0.2;        % Adaptation rate

% Parameters for other algorithms (if run)
max_iter_per_phase = 1000;
phase_num = 20;
delta_initial = 0.1;
max_iter_adaptive_sub = 20000;
max_iter_adaptive_md = 20000;
switch_step_sub = 5000;
switch_step_md = 5000;
max_iter_da = 20000;
eta_da = 0.3;
max_iter_md = 20000;
eta_md = 0.1;
max_iter_sub = 20000;
step_size_sub = 1e-3;

% Plotting Flags
plot_individual_phases = false;
plot_overall_adaptive = false; % Master plot flag for individual algorithm plots

% --- End Parameters ---

% --- Data Loading/Generation ---

load('/Users/chjiang/Dropbox/EG_EXP/Linear/WINE_Submission/Efficiency/Solver_Rating.mat', 'p_opt_solver', 'beta_opt', 'fval_solver', 'solve_time');
disp(['Loaded solver results from /Users/chjiang/Dropbox/EG_EXP/Linear/EC_Submission/Efficiency/Solver_Rating.mat']);
disp(['Solver Optimal Value (Original Func): ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
disp(['Solver Optimal Value (Original Func): ', num2str(fval_solver)]);
disp(['Solver time: ', num2str(solve_time), ' seconds']);
% p_opt_solver = p_opt_solver'; % Ensure it's a row vector

% --- Initialization ---
p_lower = max(v .* B ./ sum(abs(v),2));
p_lower(p_lower <= 0) = 1e-12;
p_upper = norm(B, 1) * ones(1, m);
mu_lower = log(p_lower);
mu_upper = log(p_upper);

L_smooth_initial = exp(max(mu_upper)) + (sum(B) / delta_initial);
sigma_calc = min(exp(mu_lower));
if sigma_calc <=0, sigma = 1e-10; else, sigma = sigma_calc; end
step_size_sub_initial = 1e-3;
gd_step_rule = @(L) 1/L;

if exist('linear_init_gd', 'file') ~= 2 || exist('linear_init_md', 'file') ~= 2
    error('Initialization functions not found.');
end
p0 = linear_init_gd(p_lower,p_upper,sum(B));
mu0 = log(p0);
x0 = linear_init_md(p0,B);

% PDHG Specific Initializations
t0_for_pdhg = ones(n, 1);
y0_for_pdhg = zeros(n, 1);

disp('--- Running Selected Algorithms ---');
results = struct();

% --- Run Fixed Step-Size PDHG ---
if run_primal_dual_pdhg_fixed
    if exist('linear_primal_dual_pdhg', 'file') == 2
        disp('Running PDHG (Fixed Step-Size)...');
        [solution_x_pdhg, solution_t_pdhg, solution_p_pdhg, solution_y_pdhg, obj_values_pdhg, dis_p_pdhg, time_pdhg, iter_pdhg] = ...
            linear_primal_dual_pdhg(v, B, x0, t0_for_pdhg, p0, y0_for_pdhg, ...
                                    max_outer_iter_pdhg, K_inner_iter_pdhg, sigma_pdhg_fixed, tau_pdhg_fixed, epsilon, ...
                                    plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['PDHG (Fixed) time: ', num2str(time_pdhg), ' seconds, Total Inner Iterations: ', num2str(iter_pdhg)]);
        results.pdhg_fixed.solution_p = solution_p_pdhg;
        results.pdhg_fixed.obj_values = obj_values_pdhg;
        results.pdhg_fixed.time = time_pdhg;
        results.pdhg_fixed.iter = iter_pdhg;
        results.pdhg_fixed.distance_final = norm(solution_p_pdhg - p_opt_solver);
    else
        warning('linear_primal_dual_pdhg.m not found. Skipping Fixed PDHG.');
        run_primal_dual_pdhg_fixed = false;
    end
end

% --- Run Adaptive Step-Size PDHG ---
if run_primal_dual_pdhg_adaptive
    if exist('linear_primal_dual_pdhg_adaptive', 'file') == 2
        disp('Running PDHG (Adaptive Step-Size)...');
        [solution_x_pdhg_a, solution_t_pdhg_a, solution_p_pdhg_a, solution_y_pdhg_a, obj_values_pdhg_a, dis_p_pdhg_a, time_pdhg_a, iter_pdhg_a] = ...
            linear_primal_dual_pdhg_adaptive(v, B, x0, t0_for_pdhg, p0, y0_for_pdhg, ...
                                             max_outer_iter_pdhg, K_inner_iter_pdhg, ...
                                             eta_initial_pdhg_adaptive, omega_initial_pdhg_adaptive, theta_pdhg_adaptive, ...
                                             epsilon, plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['PDHG (Adaptive) time: ', num2str(time_pdhg_a), ' seconds, Total Inner Iterations: ', num2str(iter_pdhg_a)]);
        results.pdhg_adaptive.solution_p = solution_p_pdhg_a;
        results.pdhg_adaptive.obj_values = obj_values_pdhg_a;
        results.pdhg_adaptive.time = time_pdhg_a;
        results.pdhg_adaptive.iter = iter_pdhg_a;
        results.pdhg_adaptive.distance_final = norm(solution_p_pdhg_a - p_opt_solver);
    else
        warning('linear_primal_dual_pdhg_adaptive.m not found. Skipping Adaptive PDHG.');
        run_primal_dual_pdhg_adaptive = false;
    end
end

% --- Run Adaptive Vanilla GD ---
if run_adaptive_gd
    if exist('linear_dual_adaptive_vanilla', 'file') == 2 && exist('linear_dual_gd', 'file') == 2
        disp('Running Adaptive Vanilla GD...');
        [solution_adpt_gd, time_adpt_gd, iter_adpt_gd, obj_values_adpt_gd, dis_adpt_gd] = ...
            linear_dual_adaptive_vanilla(v, B, mu0, max_iter_per_phase, L_smooth_initial, epsilon, mu_lower, mu_upper, delta_initial, ...
                                         plot_individual_phases, plot_overall_adaptive, p_opt_solver, fval_solver, phase_num, gd_step_rule, true);
        disp(['Adaptive Vanilla GD time: ', num2str(time_adpt_gd), ' seconds, Total Iterations: ', num2str(iter_adpt_gd)]);
        results.adaptive_gd.solution = solution_adpt_gd;
        results.adaptive_gd.obj_values = obj_values_adpt_gd;
        results.adaptive_gd.time = time_adpt_gd;
        results.adaptive_gd.iter = iter_adpt_gd;
        results.adaptive_gd.distance_final = norm(exp(solution_adpt_gd) - p_opt_solver);
    else
        warning('Required files for Adaptive GD not found. Skipping.');
        run_adaptive_gd = false;
    end
end

% --- Run APM (Adaptive NAG) ---
if run_adaptive_nag
    if exist('linear_dual_adaptive', 'file') == 2 && exist('linear_dual_agd', 'file') == 2
        disp('Running Adaptive NAG (APM)...');
        [solution_adpt_nag, time_adpt_nag, iter_adpt_nag, obj_values_adpt_nag, dis_adpt_nag] = ...
            linear_dual_adaptive(v, B, mu0, max_iter_per_phase, L_smooth_initial, sigma, epsilon, mu_lower, mu_upper, delta_initial, ...
                                 plot_individual_phases, plot_overall_adaptive, false, p_opt_solver, fval_solver, true, phase_num);
        disp(['Adaptive NAG time: ', num2str(time_adpt_nag), ' seconds, Total Iterations: ', num2str(iter_adpt_nag)]);
        results.adaptive_nag.solution = solution_adpt_nag;
        results.adaptive_nag.obj_values = obj_values_adpt_nag;
        results.adaptive_nag.time = time_adpt_nag;
        results.adaptive_nag.iter = iter_adpt_nag;
        results.adaptive_nag.distance_final = norm(exp(solution_adpt_nag) - p_opt_solver);
    else
        warning('Required files for Adaptive NAG not found. Skipping.');
        run_adaptive_nag = false;
    end
end

% --- Run Adaptive Tatonnement ---
if run_adaptive_sub
    if exist('linear_dual_subgradient_adaptive', 'file') == 2
        disp('Running Adaptive Tatonnement...');
        [solution_adpt_sub, obj_values_adpt_sub, dis_adpt_sub, time_adpt_sub, iter_adpt_sub] = ...
            linear_dual_subgradient_adaptive(v, B, p0, max_iter_adaptive_sub, step_size_sub_initial, epsilon, ...
                                             plot_overall_adaptive, p_opt_solver, fval_solver, switch_step_sub);
        disp(['Adaptive Tatonnement time: ', num2str(time_adpt_sub), ' seconds, Iterations: ', num2str(iter_adpt_sub)]);
        results.adaptive_sub.solution = solution_adpt_sub;
        results.adaptive_sub.obj_values = obj_values_adpt_sub;
        results.adaptive_sub.time = time_adpt_sub;
        results.adaptive_sub.iter = iter_adpt_sub;
        results.adaptive_sub.distance_final = norm(solution_adpt_sub - p_opt_solver);
    else
        warning('Required file for Adaptive Tatonnement not found. Skipping.');
        run_adaptive_sub = false;
    end
end

% --- Run Adaptive Mirror Descent ---
if run_adaptive_md
     if exist('linear_primal_md_adaptive', 'file') == 2
        disp('Running Adaptive Mirror Descent...');
        [solution_adpt_md, time_adpt_md, iter_adpt_md, obj_values_adpt_md, distance_adpt_md] = ...
            linear_primal_md_adaptive(v, B, x0, epsilon, max_iter_adaptive_md, ...
                                      plot_overall_adaptive, p_opt_solver, fval_solver, switch_step_md);
        disp(['Adaptive MD time: ', num2str(time_adpt_md), ' seconds, Iterations: ', num2str(iter_adpt_md)]);
        results.adaptive_md.solution = solution_adpt_md;
        results.adaptive_md.obj_values = obj_values_adpt_md;
        results.adaptive_md.time = time_adpt_md;
        results.adaptive_md.iter = iter_adpt_md;
        results.adaptive_md.distance_final = norm(solution_adpt_md - p_opt_solver);
     else
        warning('Required file for Adaptive MD not found. Skipping.');
        run_adaptive_md = false;
     end
end

% --- Run Dual Averaging ---
if run_dual_averaging
     if exist('linear_primal_da', 'file') == 2
        disp('Running Dual Averaging...');
        [solution_da, time_da, iter_da, obj_values_da, distance_da] = ...
            linear_primal_da(v, B, x0, eta_da, epsilon, max_iter_da, ...
                             plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['DA time: ', num2str(time_da), ' seconds, Iterations: ', num2str(iter_da)]);
        results.da.solution = solution_da;
        results.da.obj_values = obj_values_da;
        results.da.time = time_da;
        results.da.iter = iter_da;
        results.da.distance_final = norm(solution_da - p_opt_solver);
     else
        warning('Required file for Dual Averaging not found. Skipping.');
        run_dual_averaging = false;
     end
end

% --- Run Mirror Descent ---
if run_primal_md
    if exist('linear_primal_md', 'file') == 2
        disp('Running Mirror Descent...');
        [solution_md, time_md, iter_md, obj_values_md, distance_md] = ...
            linear_primal_md(v, B, x0, eta_md, epsilon, max_iter_md, ...
                             plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['Mirror Descent time: ', num2str(time_md), ' seconds, Iterations: ', num2str(iter_md)]);
        results.primal_md.solution = solution_md;
        results.primal_md.obj_values = obj_values_md;
        results.primal_md.time = time_md;
        results.primal_md.iter = iter_md;
        results.primal_md.distance_final = norm(solution_md - p_opt_solver);
    else
        warning('Required file for Mirror Descent not found. Skipping.');
        run_primal_md = false;
    end
end

% --- Run Tatonnement (Subgradient) ---
if run_subgradient
    if exist('linear_dual_subgradient', 'file') == 2
        disp('Running Tatonnement...');
        [solution_sub, obj_values_sub, dis_sub, time_sub, iter_sub] = ...
            linear_dual_subgradient(v, B, p0, max_iter_sub, step_size_sub, epsilon, ...
                                    plot_overall_adaptive, p_opt_solver, fval_solver);
        disp(['Tatonnement time: ', num2str(time_sub), ' seconds, Iterations: ', num2str(iter_sub)]);
        results.subgradient.solution = solution_sub;
        results.subgradient.obj_values = obj_values_sub;
        results.subgradient.time = time_sub;
        results.subgradient.iter = iter_sub;
        results.subgradient.distance_final = norm(solution_sub - p_opt_solver);
    else
        warning('Required file for Tatonnement not found. Skipping.');
        run_subgradient = false;
    end
end

disp('--- Algorithm Execution Finished ---');

% --- 绘制组合下降图 (使用原始目标差距进行比较) ---
if any([run_primal_dual_pdhg_fixed, run_primal_dual_pdhg_adaptive, run_adaptive_gd, run_adaptive_nag, run_adaptive_sub, run_adaptive_md, run_dual_averaging, run_primal_md, run_subgradient])
    disp('Plotting Objective Gap results...');
    figure;
    hold on;
    legend_entries = {};
    all_obj_values_to_plot = [];
    max_overall_iters = 0;

    % 为了一致的x轴，收集所有运行算法的最大迭代次数
    if run_adaptive_gd && isfield(results, 'adaptive_gd'), max_overall_iters = max(max_overall_iters, results.adaptive_gd.iter); end
    if run_adaptive_nag && isfield(results, 'adaptive_nag'), max_overall_iters = max(max_overall_iters, results.adaptive_nag.iter); end
    if run_adaptive_sub && isfield(results, 'adaptive_sub'), max_overall_iters = max(max_overall_iters, results.adaptive_sub.iter); end
    if run_adaptive_md && isfield(results, 'adaptive_md'), max_overall_iters = max(max_overall_iters, results.adaptive_md.iter); end
    if run_dual_averaging && isfield(results, 'da'), max_overall_iters = max(max_overall_iters, results.da.iter); end
    if run_primal_md && isfield(results, 'primal_md'), max_overall_iters = max(max_overall_iters, results.primal_md.iter); end
    if run_subgradient && isfield(results, 'subgradient'), max_overall_iters = max(max_overall_iters, results.subgradient.iter); end
    if run_primal_dual_pdhg_fixed && isfield(results, 'pdhg_fixed'), max_overall_iters = max(max_overall_iters, results.pdhg_fixed.iter); end
    if run_primal_dual_pdhg_adaptive && isfield(results, 'pdhg_adaptive'), max_overall_iters = max(max_overall_iters, results.pdhg_adaptive.iter); end
    
    % --- 绘图部分 ---

    % 绘制 Fixed PDHG 结果
    if run_primal_dual_pdhg_fixed && isfield(results, 'pdhg_fixed') && ~isempty(results.pdhg_fixed.obj_values)
        total_iters = results.pdhg_fixed.iter;
        semilogy(0:total_iters-1, abs(results.pdhg_fixed.obj_values), '-<', 'DisplayName', 'PDHG Fixed', 'LineWidth', 2);
        legend_entries{end+1} = 'PDHG Fixed';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.pdhg_fixed.obj_values(:))];
    end

    % 绘制 Adaptive PDHG 结果
    if run_primal_dual_pdhg_adaptive && isfield(results, 'pdhg_adaptive') && ~isempty(results.pdhg_adaptive.obj_values)
        total_iters = results.pdhg_adaptive.iter;
        semilogy(0:total_iters-1, abs(results.pdhg_adaptive.obj_values), '->', 'DisplayName', 'PDHG Adaptive', 'LineWidth', 2);
        legend_entries{end+1} = 'PDHG Adaptive';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.pdhg_adaptive.obj_values(:))];
    end

    % 绘制 Adaptive GD 结果
    if run_adaptive_gd && isfield(results, 'adaptive_gd') && ~isempty(results.adaptive_gd.obj_values)
        total_iters = results.adaptive_gd.iter;
        semilogy(0:total_iters-1, abs(results.adaptive_gd.obj_values), '-^', 'DisplayName', 'Adaptive GD', 'LineWidth', 2);
        legend_entries{end+1} = 'Adaptive GD';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_gd.obj_values(:))];
    end

    % 绘制 Tatonnement 结果
    if run_subgradient && isfield(results, 'subgradient') && ~isempty(results.subgradient.obj_values)
        total_iters = results.subgradient.iter;
        semilogy(0:total_iters-1, abs(results.subgradient.obj_values), '-v', 'DisplayName', 'Tatonnement', 'LineWidth', 2);
        legend_entries{end+1} = 'Tatonnement';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.subgradient.obj_values(:))];
    end

    % 绘制 Adaptive Tatonnement 结果
    if run_adaptive_sub && isfield(results, 'adaptive_sub') && ~isempty(results.adaptive_sub.obj_values)
        total_iters = results.adaptive_sub.iter;
        semilogy(0:total_iters-1, abs(results.adaptive_sub.obj_values), '-s', 'DisplayName', 'Adaptive Tatonnement', 'LineWidth', 2);
        legend_entries{end+1} = 'Adaptive Tatonnement';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_sub.obj_values(:))];
    end

    % 绘制 Mirror Descent 结果
    if run_primal_md && isfield(results, 'primal_md') && ~isempty(results.primal_md.obj_values)
        total_iters = results.primal_md.iter;
        semilogy(0:total_iters-2, abs(results.primal_md.obj_values), '-x', 'DisplayName', 'Mirror Descent', 'LineWidth', 2);
        legend_entries{end+1} = 'Mirror Descent';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.primal_md.obj_values(:))];
    end


    % 绘制 Adaptive MD 结果
    if run_adaptive_md && isfield(results, 'adaptive_md') && ~isempty(results.adaptive_md.obj_values)
        total_iters = results.adaptive_md.iter;
        semilogy(0:total_iters-2, abs(results.adaptive_md.obj_values), '-o', 'DisplayName', 'Adaptive MD', 'LineWidth', 2);
        legend_entries{end+1} = 'Adaptive MD';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_md.obj_values(:))];
    end


     % 绘制 Dual Averaging 结果
     if run_dual_averaging && isfield(results, 'da') && ~isempty(results.da.obj_values)
        total_iters = results.da.iter;
        semilogy(0:total_iters, abs(results.da.obj_values), '-*', 'DisplayName', 'Dual Averaging', 'LineWidth', 2);
        legend_entries{end+1} = 'PDHG';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.da.obj_values(:))];
    end

    % 绘制 Adaptive NAG 结果
    if run_adaptive_nag && isfield(results, 'adaptive_nag') && ~isempty(results.adaptive_nag.obj_values)
        total_iters = results.adaptive_nag.iter;
        semilogy(0:total_iters-1, abs(results.adaptive_nag.obj_values), '-d', 'DisplayName', 'APM', 'LineWidth', 2);
        legend_entries{end+1} = 'APM';
        all_obj_values_to_plot = [all_obj_values_to_plot; abs(results.adaptive_nag.obj_values(:))];
    end
    
    hold off;

    % 自定义绘图
    % Customize plot
    set(gca, 'FontSize', 15); % Axis font size
    xlabel('Iteration', 'FontSize', 20);
    ylabel('Objective Value Gap', 'FontSize', 20);
    title_str_obj = sprintf('n=%d, m=%d', n, m);
    title(title_str_obj, 'FontSize', 25);
    
    % Set logarithmic y-axis with proper ticks
    set(gca, 'YScale', 'log');
    grid on;
    
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
        
        % Create logarithmic ticks
        min_exp = floor(log10(min_val));
        max_exp = ceil(log10(max_val));
        yticks(10.^(min_exp:1:max_exp));
        yticklabels(arrayfun(@(x) sprintf('10^{%d}', x), min_exp:1:max_exp, 'UniformOutput', false));
    else
        ylim([epsilon/10, 1]);
        yticks(10.^(-2:1:0));
        yticklabels({'10^{-2}', '10^{-1}', '10^{0}'});
    end

    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'northeast', 'FontSize', 10);
    else
        warning('No algorithms produced results for objective gap plotting.');
    end
else
    disp('No algorithms produced results for objective gap plotting.');
end

% --- 与最优价格 p* 的最终距离 ---
disp('--- Final Distance to Optimal Primal Solution p* ---');
found_results = false;
if run_primal_dual_pdhg_fixed && isfield(results, 'pdhg_fixed')
    disp(['PDHG (Fixed): ', num2str(results.pdhg_fixed.distance_final)]);
    found_results = true;
end
if run_primal_dual_pdhg_adaptive && isfield(results, 'pdhg_adaptive')
    disp(['PDHG (Adaptive): ', num2str(results.pdhg_adaptive.distance_final)]);
    found_results = true;
end
if run_adaptive_gd && isfield(results, 'adaptive_gd'), disp(['Adaptive GD: ', num2str(results.adaptive_gd.distance_final)]); found_results = true; end
if run_adaptive_nag && isfield(results, 'adaptive_nag'), disp(['APM: ', num2str(results.adaptive_nag.distance_final)]); found_results = true; end
if run_adaptive_sub && isfield(results, 'adaptive_sub'), disp(['Adaptive Tatonnement: ', num2str(results.adaptive_sub.distance_final)]); found_results = true; end
if run_adaptive_md && isfield(results, 'adaptive_md'), disp(['Adaptive MD: ', num2str(results.adaptive_md.distance_final)]); found_results = true; end
if run_dual_averaging && isfield(results, 'da'), disp(['Dual Averaging: ', num2str(results.da.distance_final)]); found_results = true; end
if run_primal_md && isfield(results, 'primal_md'), disp(['Mirror Descent: ', num2str(results.primal_md.distance_final)]); found_results = true; end
if run_subgradient && isfield(results, 'subgradient'), disp(['Tatonnement: ', num2str(results.subgradient.distance_final)]); found_results = true; end
if ~found_results
    disp('No results to display.');
end

disp('--- Script Finished ---');