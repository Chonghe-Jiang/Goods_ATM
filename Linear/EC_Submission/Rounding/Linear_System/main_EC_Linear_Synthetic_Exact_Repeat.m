clc
clear
addpath('/Users/chjiang/Dropbox/EG_EXP/Linear/EC_Submission/Efficiency');

%%% Define the sizes to test
sizes = [50, 100, 200, 300, 400];
num_repeats = 10;  % Number of times to repeat the experiment for each size

%%% Define the CSV filename
csv_filename = 'solver_adaptive_gap_data_repeat.csv';

%%% Loop through each size
for size_idx = 1:length(sizes)
    n = sizes(size_idx);  % Number of rows
    m = sizes(size_idx);  % Number of columns
    
    % Loop through each repeat
    for repeat = 1:num_repeats
        % Generate random B vector
        B = ones(n, 1);
        
        % Generate random v matrix
        v = randi([1, 10], n, m);
        
        % Define solver parameters
        max_iter = 10000;
        max_iter_adaptive = 4500;
        p_lower = max(v .* B ./ sum(abs(v), 2));
        p_upper = norm(B, 1) * ones(1, m);
        mu_lower = log(p_lower);
        mu_upper = log(p_upper);
        delta = 0.1;
        epsilon = 0.1;
        sigma = min(exp(mu_lower));
        L = exp(max(mu_upper)) + (sum(B) / delta);
        p0 = linear_init_gd(p_lower, p_upper, sum(B));
        mu0 = log(p0);
        
        % Solve the problem using the solver
        tic;
        [p_opt_solver, beta_opt, fval_solver, solve_time] = linear_dual_solver(n, m, B, v);
        solver_time = toc;
        
        % Calculate the Gap
        gap_optimal = linear_optimal_gap(v, p_opt_solver');
        
        % Run the adaptive algorithm
        adaptive_plot_flag = false;
        plot_flag = false;
        plot_flag_smooth = false;
        adaptive = true;
        phase_num = 200;
        tic;
        [solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_values_adaptive, dis_adaptive, results_matrix] = linear_dual_adaptive_exact(v, B, mu0, max_iter_adaptive, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver', fval_solver, adaptive, phase_num);
        adaptive_time = toc;
        
        % Calculate the solution gap between Exact and Solver
        Exact_solver_gap = norm(solution_adaptive - log(p_opt_solver'));
        
        % Create a table with the data, including total_iter_adaptive as the last column
        data_table = table(n, m, solver_time, adaptive_time, gap_optimal, Exact_solver_gap, total_iter_adaptive, ...
            'VariableNames', {'n', 'm', 'solver_time', 'adaptive_time', 'gap_optimal', 'Exact_solver_gap', 'total_iter_adaptive'});
        
        % Append the data to the CSV file
        if exist(csv_filename, 'file') == 2
            % File exists, append the data
            writetable(data_table, csv_filename, 'WriteMode', 'append');
        else
            % File does not exist, create a new file with headers
            writetable(data_table, csv_filename);
        end
        
        % Display the results
        disp(['Size: ', num2str(n), 'x', num2str(m), ', Repeat: ', num2str(repeat)]);
        disp(['Solver time: ', num2str(solver_time), ' seconds']);
        disp(['Adaptive AGD time: ', num2str(adaptive_time), ' seconds']);
        disp(['Gap: ', num2str(gap_optimal)]);
        disp(['Exact Solver Gap: ', num2str(Exact_solver_gap)]);
        disp(['Total Adaptive Iterations: ', num2str(total_iter_adaptive)]);
        disp('----------------------------------------');
    end
end

disp(['All data saved to ', csv_filename]);