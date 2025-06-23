function [solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_adaptive_vanilla, dis_adaptive_vanilla] = linear_dual_adaptive_vanilla(v, B, mu_0, max_iter_per_phase, L_initial, epsilon, mu_lower, mu_upper, delta_initial, plot_flag, adaptive_plot_flag, p_opt_solver, fval_solver, phase_num, step_size_rule, adaptive)
    % Adaptive strategy using Vanilla Gradient Descent (GD) in phases.
    % Input:
    % v - parameter matrix v \in R^{n*m}
    % B - vector B \in R^{n*1}
    % mu_0 - initial point mu_0 \in R^{1*m}
    % max_iter_per_phase - maximum number of iterations *per phase*
    % L_initial - initial Lipschitz constant estimate (for smoothed gradient)
    % epsilon - tolerance for stopping criteria (based on original objective value gap)
    % mu_lower - lower bound for mu
    % mu_upper - upper bound for mu
    % delta_initial - initial smoothing parameter
    % plot_flag - flag to plot results for each individual GD phase call (usually false)
    % adaptive_plot_flag - flag to plot overall adaptive results
    % p_opt_solver - optimal solution (primal) from the solver
    % fval_solver - optimal function value (original function) from the solver
    % phase_num - total number of phases
    % step_size_rule - Function handle or value for GD step size (e.g., @(L) 1/L or a fixed value)

    % Output:
    % solution_adaptive - final solution (mu) after all phases
    % total_time_adaptive - total time taken across all phases
    % total_iter_adaptive - total iterations accumulated across all phases
    % obj_adaptive_vanilla - array of original objective function value gaps across all iterations
    % dis_adaptive_vanilla - array of distances to p_opt_solver across all iterations

    % Initialize tracking variables
    obj_adaptive_vanilla = [];
    dis_adaptive_vanilla = [];
    total_time_adaptive = 0;
    total_iter_adaptive = 0;

    % Initialize phase-dependent parameters
    mu_current = mu_0;
    delta = delta_initial;
    L = L_initial; % L corresponds to the smoothed function's gradient Lipschitz constant

    fprintf('--- Starting Adaptive Vanilla GD ---\n');
    % Loop over phases
    for phase = 1:phase_num
        fprintf('--- Phase %d/%d ---\n', phase, phase_num);
        fprintf('Delta = %e, L = %e\n', delta, L);

        % Determine step size for this phase based on the rule
        if isa(step_size_rule, 'function_handle')
            step_size = step_size_rule(L);
        else
            step_size = step_size_rule; % Use fixed value if provided
        end
        fprintf('Step Size = %e\n', step_size);

        % Call the linear_dual_gd function for the current phase
        % Note: linear_dual_gd returns obj_values_gd (original gap), f_smooth_values_gd, dis_gd
        [solution_phase, time_phase, iter_phase, obj_phase_gd, ~, dis_phase_gd, convergence] = ...
            linear_dual_gd(v, B, mu_current, max_iter_per_phase, L, epsilon, mu_lower, mu_upper, delta, ...
                           plot_flag, p_opt_solver, fval_solver, step_size, adaptive);

        % Append results from this phase
        obj_adaptive_vanilla = [obj_adaptive_vanilla; obj_phase_gd];
        dis_adaptive_vanilla = [dis_adaptive_vanilla; dis_phase_gd];

        % Accumulate total time and iterations
        total_time_adaptive = total_time_adaptive + time_phase;
        total_iter_adaptive = total_iter_adaptive + iter_phase;

        % Print phase summary
        fprintf('Phase %d: Iterations = %d, Time = %.4f s, Final Orig. Obj. Gap = %e\n', ...
                phase, iter_phase, time_phase, obj_phase_gd(end));

        % Update starting point for the next phase
        mu_current = solution_phase;
        solution_adaptive = mu_current; % Store the latest solution

        % Check for overall convergence (based on the last phase's result)
        % Using the original objective gap check from linear_dual_gd
        if convergence || obj_phase_gd(end) < epsilon
            fprintf('Convergence tolerance reached in Phase %d.\n', phase);
            break; % Exit the phase loop
        end

        % If not the last phase, update delta and L for the next phase
        if phase < phase_num
            % Update delta and L using the same logic as linear_dual_adaptive
            % You might want to adjust this strategy for vanilla GD
            if phase >=1 && phase <=3
                delta = delta / 3;
            elseif phase <= 5
                delta = delta / 3; % Original code had this twice, kept for consistency
            else
                delta = delta / 1.5;
            end
            % Recalculate L based on the new delta
            L = exp(max(mu_upper)) + sum(B) / delta;
        end
    end % End of phase loop

    fprintf('--- Adaptive Vanilla GD Finished ---\n');
    fprintf('Total Iterations = %d, Total Time = %.4f s\n', total_iter_adaptive, total_time_adaptive);

    % Plot the overall results if adaptive_plot_flag is true
    if adaptive_plot_flag && ~isempty(obj_adaptive_vanilla)
        figure;
        iterations_total = 1:total_iter_adaptive; % Correct x-axis for accumulated iterations

        subplot(2, 1, 1);
        semilogy(iterations_total, abs(obj_adaptive_vanilla), 'LineWidth', 1.5);
        xlabel('Total Iterations');
        ylabel('Absolute Original Objective Gap');
        title('Adaptive Vanilla GD - Objective Value Convergence');
        grid on;
        if convergence || (~isempty(obj_adaptive_vanilla) && obj_adaptive_vanilla(end) < epsilon)
             hold on;
             semilogy(total_iter_adaptive, abs(obj_adaptive_vanilla(end)), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Mark convergence point
             hold off;
        end


        subplot(2, 1, 2);
        plot(iterations_total, dis_adaptive_vanilla, 'LineWidth', 1.5);
        xlabel('Total Iterations');
        ylabel('Distance ||exp(mu) - p*||');
        title('Adaptive Vanilla GD - Distance to Optimal Solution');
        grid on;
         if convergence || (~isempty(obj_adaptive_vanilla) && obj_adaptive_vanilla(end) < epsilon)
             hold on;
             plot(total_iter_adaptive, dis_adaptive_vanilla(end), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Mark convergence point
             hold off;
        end
    end
end