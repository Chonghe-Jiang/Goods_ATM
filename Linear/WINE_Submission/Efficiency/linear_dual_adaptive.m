function [solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_adaptive, dis_adaptive] = linear_dual_adaptive(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth,  p_opt_solver, fval_solver, adaptive, phase_num)
    % Adaptive NAG (Nesterov's Accelerated Gradient) strategy with multiple phases.
    % Input:
    % v                 - parameter matrix v \in R^{n*m}
    % B                 - vector B \in R^{n*1}
    % mu_0              - initial point mu_0 \in R^{1*m}
    % max_iter          - maximum number of iterations *per phase*
    % L                 - initial Lipschitz constant estimate
    % sigma             - initial strong convexity parameter estimate
    % epsilon           - tolerance for stopping criteria (based on original objective value gap)
    % mu_lower          - lower bound for mu
    % mu_upper          - upper bound for mu
    % delta             - initial smoothing parameter
    % plot_flag         - flag to plot results for each individual AGD phase call
    % adaptive_plot_flag- flag to plot overall adaptive results from this function
    % plot_flag_smooth  - flag to plot smoothed objective separately from within AGD calls
    % p_opt_solver      - optimal solution (primal) from the solver
    % fval_solver       - optimal function value (original function) from the solver
    % adaptive          - boolean flag passed to AGD to enable its internal adaptive checks
    % phase_num         - total number of phases

    % Output:
    % solution_adaptive   - final solution (mu) after all phases
    % total_time_adaptive - total time taken across all phases
    % total_iter_adaptive - total iterations accumulated across all phases
    % obj_adaptive        - array of original objective function value gaps across all iterations
    % dis_adaptive        - array of distances to p_opt_solver across all iterations

    % --- Initialization ---
    % Initialize output variables *before* the loop
    solution_adaptive = mu_0; % Initialize with the starting point
    total_time_adaptive = 0;
    total_iter_adaptive = 0;
    obj_adaptive = []; % Will be concatenated inside the loop
    dis_adaptive = []; % Will be concatenated inside the loop

    % Initialize phase-dependent parameters
    mu_current = mu_0; % Use a separate variable for the current phase's start
    delta_current = delta;
    L_current = L;
    % --- End Initialization ---

    % Optional: Check if phase_num is valid
    if phase_num < 1
        warning('linear_dual_adaptive:NoPhases', 'phase_num is less than 1. No phases will be run. Returning initial state.');
        % Ensure outputs that should be arrays are empty if no loop runs
        obj_adaptive = [];
        dis_adaptive = [];
        return; % Exit the function early
    end

    fprintf('--- Starting Adaptive NAG ---\n');
    % --- Loop over phases ---
    for phase = 1:phase_num
        fprintf('--- Phase %d/%d ---\n', phase, phase_num);
        fprintf('Delta = %e, L = %e\n', delta_current, L_current);

        % Call the linear_dual_agd function for the current phase
        % Make sure linear_dual_agd returns f_smooth_values_agd if needed elsewhere
        [solution_phase, time_phase, iter_phase, obj_phase, ~, dis_phase, convergence] = ...
            linear_dual_agd(v, B, mu_current, max_iter, L_current, sigma, epsilon, mu_lower, mu_upper, delta_current, ...
                            plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive);

        % Append results from this phase
        obj_adaptive = [obj_adaptive; obj_phase];
        dis_adaptive = [dis_adaptive; dis_phase];

        % Accumulate total time and iterations
        total_time_adaptive = total_time_adaptive + time_phase;
        total_iter_adaptive = total_iter_adaptive + iter_phase;

        % Print phase summary
        fprintf('Phase %d: Iterations = %d, Time = %.4f s, Final Orig. Obj. Gap = %e\n', ...
                phase, iter_phase, time_phase, obj_phase(end));

        % Update starting point for the next phase
        mu_current = solution_phase;
        solution_adaptive = mu_current; % Store the latest valid solution

        % Check for overall convergence (based on the primary criterion)
        if convergence % Check the flag returned by linear_dual_agd
            fprintf('Convergence tolerance (epsilon) reached in Phase %d.\n', phase);
            break; % Exit the phase loop
        end

        % If not the last phase, update delta and L for the next phase
        if phase < phase_num
            % Update delta and L using the defined strategy
            if phase >=1 && phase <=3
                delta_current = delta_current / 5;
            elseif phase <= 5
                delta_current = delta_current / 3; % Original code had this twice
            else
                delta_current = delta_current / 1.5;
            end
            % Recalculate L based on the new delta
            L_current = exp(max(mu_upper)) + sum(B) / delta_current;
            % Note: Sigma is not updated in this scheme, might need adjustment if desired
        end
    end % End of phase loop
    % --- End Loop ---

    fprintf('--- Adaptive NAG Finished ---\n');
    fprintf('Total Iterations = %d, Total Time = %.4f s\n', total_iter_adaptive, total_time_adaptive);

    % --- Plotting Overall Results ---
    if adaptive_plot_flag && ~isempty(obj_adaptive)
        figure; % Create a new figure for the overall adaptive results
        iterations_total = 1:total_iter_adaptive; % Correct x-axis for accumulated iterations

        subplot(2, 1, 1);
        semilogy(iterations_total, abs(obj_adaptive), 'LineWidth', 1.5);
        xlabel('Total Iterations');
        ylabel('Original Objective Gap |F(\cdot) - F*|');
        title('Adaptive NAG - Overall Objective Value Convergence');
        grid on;
        hold on; yline(epsilon, 'r--', 'Tolerance'); hold off; % Show tolerance

        subplot(2, 1, 2);
        plot(iterations_total, dis_adaptive, 'LineWidth', 1.5);
        xlabel('Total Iterations');
        ylabel('Distance ||exp(mu) - p*||');
        title('Adaptive NAG - Overall Distance to Optimal Solution');
        grid on;
    end
    % --- End Plotting ---
end
