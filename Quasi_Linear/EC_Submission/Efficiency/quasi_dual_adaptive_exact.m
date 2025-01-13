function [solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_adaptive, dis_adaptive, results_matrix] = quasi_dual_adaptive_exact(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive, phase_num)
    % Initialize variables
    time_adaptive = [];
    iter_adaptive = [];
    obj_adaptive = [];
    dis_adaptive = [];
    total_time_adaptive = 0;  % Initialize total time
    total_iter_adaptive = 0;
    results_matrix = [];  % Initialize the results matrix

    % Loop over phases
    for phase = 1:phase_num
        % Call the linear_dual_agd function
        %%% ! Step 1: Run the inner step
        [solution_phase, time_phase, iter_phase, obj_phase, dis_phase] = quasi_dual_agd_exact(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive);
        %%% ! Step 2: Extract the gap_current as the paper mentioned - functino value not variable value
        % Calculate and print the gap between solution_phase and log(p_opt_solver)
        % gap_current = norm(solution_phase - log(p_opt_solver), 2);
        % fprintf('Phase %d: Gap between solution and log(p_opt_solver) = %f\n', phase, gap_current);
        
        % Update variables
        time_adaptive = [time_adaptive; time_phase];
        iter_adaptive = [iter_adaptive; iter_phase];
        obj_adaptive = [obj_adaptive; obj_phase];
        dis_adaptive = [dis_adaptive; dis_phase];
        
        % Accumulate total time
        total_time_adaptive = total_time_adaptive + time_phase;
        total_iter_adaptive = total_iter_adaptive + iter_phase;
        %%% ! Print the information
        % Print phase information
        fprintf('Phase %d: Iterations = %d, Initial Objective = %f, Final Objective = %f\n', phase, iter_phase, obj_phase(1), obj_phase(end));
        gap_current = obj_phase(end);
        %%% ! Step3: Set the radius according to the paper
        % Call the exact oracle function
        radius_curent = sqrt(2*gap_current/sigma);
        %%% ! Step4: Run the oracle
        [oracle_result, optimal_ce, optimal_value, sum_exp_mu] = quasi_exact_oracle(v, B, solution_phase, radius_curent);
        
        % Store the results in the matrix
        results_matrix = [results_matrix; optimal_value, sum_exp_mu];
        
        % Check for convergence using the oracle result
        if oracle_result
            fprintf('Oracle result is True. Convergence achieved.\n');
            solution_adaptive = optimal_ce;
            break;
        end
        
        %%% Todo-Check whether this is the most efficient way
        if phase >= 1 && phase <= 3
            delta = delta / 3;
            L = exp(max(mu_upper)) + sum(B) / delta;
            mu_0 = solution_phase;
        elseif phase <= 5
            delta = delta / 2;
            L = exp(max(mu_upper)) + sum(B) / delta;
            mu_0 = solution_phase;
        else
            delta = delta / 1.5;
            L = exp(max(mu_upper)) + sum(B) / delta;
            mu_0 = solution_phase;
        end
    end

    % If convergence was not achieved, set the final solution to the last phase solution
    if ~exist('solution_adaptive', 'var')
        solution_adaptive = solution_phase;
    end

    % Plot the results if plot_flag is true
    if adaptive_plot_flag
        figure;
        iterations = 1:length(obj_adaptive);
        
        subplot(2, 1, 1);
        plot(iterations, abs(obj_adaptive), 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Function Value Gap');
        title('Adaptive SAG - Function Value Convergence');
        grid on;
        
        subplot(2, 1, 2);
        plot(iterations, dis_adaptive, 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Iteration Distance');
        title('Adaptive SAG - Iteration Convergence');
        grid on;
    end
end