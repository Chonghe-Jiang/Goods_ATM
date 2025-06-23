function [solution, time, iter, obj_values_gd, f_smooth_values_gd, dis_gd, convergence] = linear_dual_gd(v, B, mu_0, max_iter, L, epsilon, mu_lower, mu_upper, delta, plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, step_size, adaptive)
    % Solves the dual problem using vanilla Gradient Descent on the smoothed objective.
    % Includes adaptive phase termination criterion (if adaptive=true).
    % Includes option to plot smoothed objective separately (plot_flag_smooth).
    % Organization aligned with linear_dual_agd for consistency.
    % Input:
    % v                - parameter matrix v \in R^{n*m}
    % B                - vector B \in R^{n*1}
    % mu_0             - initial point mu_0 \in R^{1*m}
    % max_iter         - maximum number of iterations
    % L                - Lipschitz constant of the smoothed function gradient
    % epsilon          - tolerance for stopping criteria (based on original objective value gap)
    % mu_lower         - lower bound for mu
    % mu_upper         - upper bound for mu
    % delta            - smoothing parameter
    % plot_flag        - flag to plot primary results (original gap, smooth value, distance)
    % plot_flag_smooth - flag to plot smoothed objective value in a separate figure
    % p_opt_solver     - optimal solution (primal) from the mosek solver
    % fval_solver      - optimal function value (original function) from the mosek solver
    % step_size        - step size for gradient descent (e.g., 1/L)
    % adaptive         - boolean flag to enable adaptive phase termination check

    % Output:
    % solution           - solution (mu) of the optimization problem
    % time               - time taken to solve the problem
    % iter               - number of iterations performed (actual count)
    % obj_values_gd      - array of original objective function values gap (f(mu) - fval_solver)
    % f_smooth_values_gd - array of smoothed objective function values
    % dis_gd             - array of distances between exp(mu) at each iteration and p_opt_solver
    % convergence        - boolean flag indicating if epsilon convergence was met (primary criterion)

    % --- Input Validation ---
    % Validate adaptive flag
    if nargin < 15 % Check number of arguments based on function definition
        adaptive = false; % Default to false if not provided
    end
    if ~islogical(adaptive) && ~isnumeric(adaptive)
        error('Input ''adaptive'' must be a logical value (true/false) or numeric (1/0).');
    end
    adaptive = logical(adaptive); % Ensure it's logical

    % Validate plot_flag_smooth
    if nargin < 11 || isempty(plot_flag_smooth) % Check if plot_flag_smooth was provided
        plot_flag_smooth = false; % Default to false if not provided
    end
     if ~islogical(plot_flag_smooth) && ~isnumeric(plot_flag_smooth)
         error('Input ''plot_flag_smooth'' must be a logical value (true/false) or numeric (1/0).');
    end
    plot_flag_smooth = logical(plot_flag_smooth); % Ensure it's logical

    % Validate plot_flag
     if nargin < 10 || isempty(plot_flag)
        plot_flag = false; % Default plot_flag
     end
     if ~islogical(plot_flag) && ~isnumeric(plot_flag)
         error('Input ''plot_flag'' must be a logical value (true/false) or numeric (1/0).');
     end
     plot_flag = logical(plot_flag);

    % Determine step size if not provided
    if nargin < 14 || isempty(step_size) % step_size is the 14th argument now
        step_size = 1 / L; % Default step size based on smoothed Lipschitz constant
        disp(['Using default GD step size: ', num2str(step_size)]);
    elseif step_size <= 0
        error('Step size must be positive.');
    end
    % --- End Input Validation ---

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Define the projection operator
    P = @(mu) max(mu_lower, min(mu, mu_upper));

    % Initialize the variables
    mu = mu_0;
    mu_new = mu_0; % Initialize mu_new

    % Initialize arrays to store objective function values and distances
    obj_values_gd = zeros(max_iter, 1);
    f_smooth_values_gd = zeros(max_iter, 1);
    dis_gd = zeros(max_iter, 1);
    convergence = false; % Flag for primary epsilon convergence

    % Start timing
    tic;

    final_iter = 0; % Track the actual number of iterations performed
    % --- Gradient Descent Iterations ---
    for iter = 1:max_iter
        final_iter = iter; % Update final iteration count

        % Compute the original objective function value gap
        obj_orig = sum(exp(mu)) + sum(B .* max(log(v)-mu, [], 2)) - fval_solver;

        % Compute the smoothed objective function value f_smooth(mu)
        max_log_v_mu = max(log(v) - repmat(mu, n, 1), [], 2);
        rescaled_log_v_mu = (log(v) - repmat(mu, n, 1) - max_log_v_mu) / delta;
        log_sum_exp_term = log(sum(exp(rescaled_log_v_mu), 2));
        log_sum_exp_term(isinf(log_sum_exp_term)) = -realmax; % Handle potential -Inf
        f_smooth = sum(exp(mu)) + delta * sum(B .* ((max_log_v_mu/delta) + log_sum_exp_term));

        % Store tracking values for the current iteration 'iter' based on state 'mu'
        obj_values_gd(iter) = obj_orig;
        f_smooth_values_gd(iter) = f_smooth;
        dis_gd(iter) = norm(exp(mu) - p_opt_solver);

        % Calculate the gradient of the smoothed function f_smooth at the current mu
        temp1 = log(v) - repmat(mu, n, 1); % Using mu
        max_temp1 = max(temp1,[],2);
        exp_temp1 = exp((temp1 - max_temp1)/delta);
        sum_exp_temp1 = sum(exp_temp1, 2);
        valid_rows = sum_exp_temp1 > 1e-100; % Use a small threshold for stability
        cal_temp1 = zeros(n, m);
        if any(valid_rows)
            cal_temp1(valid_rows, :) = exp_temp1(valid_rows, :) ./ sum_exp_temp1(valid_rows);
        end
        temp_2 = sum(B .* cal_temp1);
        grad_f_smooth = exp(mu) - temp_2; % Gradient of smoothed function at mu

        % Update mu using Gradient Descent step
        mu_new = P(mu - step_size * grad_f_smooth);

        % --- Stopping Criteria ---
        % 1. Primary Convergence Check (Original Objective Gap)
        if iter >= 2 && obj_orig < epsilon
             convergence = true; % Set primary convergence flag
             break; % Exit loop
        end

        % 2. Adaptive Phase Termination Check (Smoothed Objective Plateau)
        %    Only active if 'adaptive' flag is true.
        if adaptive && iter >= 200 % Requires iter >= 200 (as in AGD)
            if abs(f_smooth_values_gd(iter) - f_smooth_values_gd(iter-1)) < 1e-3 && ...
               abs(f_smooth_values_gd(iter-1) - f_smooth_values_gd(iter-2)) < 1e-3 && ...
               abs(f_smooth_values_gd(iter-2) - f_smooth_values_gd(iter-3)) < 1e-3 && ...
               abs(f_smooth_values_gd(iter-3) - f_smooth_values_gd(iter-4)) < 1e-3
                fprintf('GD Info: Adaptive phase termination triggered at iteration %d.\n', iter);
                break; % Exit loop (phase termination, not final convergence)
            end
        end
        % --- End Stopping Criteria ---

        % Update mu for the next iteration
        mu = mu_new;

        % Check for NaN/Inf in mu (optional safeguard)
        if any(isnan(mu)) || any(isinf(mu))
            warning('GD:NaNInf', 'NaN or Inf detected in mu at iteration %d. Stopping.', iter);
            break; % Exit loop
        end

    end % End of iterations loop
    % --- End Gradient Descent Iterations ---

    % End timing
    time = toc;

    % Extract the solution - use the state 'mu' which corresponds to the last calculated objectives/gradient
    solution = mu;

    % Trim the arrays to the actual number of iterations performed
    obj_values_gd = obj_values_gd(1:final_iter);
    f_smooth_values_gd = f_smooth_values_gd(1:final_iter);
    dis_gd = dis_gd(1:final_iter);

    % Final convergence check in case max_iter was reached exactly when obj_orig < epsilon
    if ~convergence && ~isempty(obj_values_gd) && obj_values_gd(end) < epsilon
        convergence = true;
    end

    if ~convergence && final_iter == max_iter
         fprintf('GD Warning: Did not converge within %d iterations. Final objective gap: %e\n', max_iter, obj_values_gd(end));
    end

    % --- Plotting ---
    % Plot the primary results if plot_flag is true
    if plot_flag
        figure; % Create figure for primary plots
        plot_iterations = 1:final_iter;

        subplot(3, 1, 1);
        semilogy(plot_iterations, abs(obj_values_gd), 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('Original Objective Gap');
        title('GD - Original Objective Value Convergence');
        grid on;
        hold on; yline(epsilon, 'r--', 'Tolerance'); hold off; % Show tolerance line

        subplot(3, 1, 2);
        plot(plot_iterations, f_smooth_values_gd, 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('Smoothed Objective Value');
        title('GD - Smoothed Objective Value Convergence');
        grid on;

        subplot(3, 1, 3);
        plot(plot_iterations, dis_gd, 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('Distance ||exp(mu) - p*||');
        title('GD - Distance to Optimal Solution (Primal Space)');
        grid on;
    end

    % Plot smoothed objective function value separately if plot_flag_smooth is true
    if plot_flag_smooth
        figure; % Create a new figure for the smooth plot
        plot_iterations = 1:final_iter;
        plot(plot_iterations, f_smooth_values_gd, 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('Smoothed Objective Value');
        title('GD - Smoothed Objective Function Value');
        grid on;
    end
    % --- End Plotting ---

end % End of function
