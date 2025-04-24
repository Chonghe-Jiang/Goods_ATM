function [solution, time, iter, obj_values, f_smooth_values_agd, dis_agd, convergence] = linear_dual_agd(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive)
    % Input:
    % v                 - parameter matrix v \in R^{n*m}
    % B                 - vector B \in R^{n*1}
    % mu_0              - initial point mu_0 \in R^{1*m}
    % max_iter          - maximum number of iterations
    % L                 - Lipschitz constant
    % sigma             - strong convexity parameter
    % epsilon           - tolerance
    % mu_lower          - lower bound for mu
    % mu_upper          - upper bound for mu
    % delta             - smoothing parameter
    % plot_flag         - flag to plot primary results (original gap, distance)
    % plot_flag_smooth  - flag to plot smoothed objective value in a separate figure
    % p_opt_solver      - optimal solution from the mosek solver
    % fval_solver       - optimal function value from the mosek solver
    % adaptive          - boolean flag to enable adaptive phase termination check

    % Output:
    % solution            - solution of the optimization problem
    % time                - time taken to solve the problem
    % iter                - actual number of iterations performed
    % obj_values          - array of original objective function values gap
    % f_smooth_values_agd - array of smoothed objective function values (NEW OUTPUT)
    % dis_agd             - array of distances between exp(mu) at each iteration and p_opt_solver
    % convergence         - boolean flag indicating if epsilon convergence was met

    % --- Input Validation ---
    if nargin < 15
        adaptive = false;
    end
    if ~islogical(adaptive) && ~isnumeric(adaptive)
        error('Input ''adaptive'' must be a logical value (true/false) or numeric (1/0).');
    end
    adaptive = logical(adaptive);

    if nargin < 12 || isempty(plot_flag_smooth)
        plot_flag_smooth = false;
    end
     if ~islogical(plot_flag_smooth) && ~isnumeric(plot_flag_smooth)
         error('Input ''plot_flag_smooth'' must be a logical value (true/false) or numeric (1/0).');
     end
    plot_flag_smooth = logical(plot_flag_smooth);

    if nargin < 11 || isempty(plot_flag)
        plot_flag = false;
    end
     if ~islogical(plot_flag) && ~isnumeric(plot_flag)
         error('Input ''plot_flag'' must be a logical value (true/false) or numeric (1/0).');
     end
    plot_flag = logical(plot_flag);
    % --- End Input Validation ---

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Define the projection operator
    P = @(mu) max(mu_lower, min(mu, mu_upper));

    % Initialize the variables and parameters
    mu = mu_0;
    y = mu_0;
    mu_new = mu_0; % Initialize mu_new
    y_new = y;     % Initialize y_new
    % q = 0.1;
    %%% Todo - choose the best parameter here
    if L <= 0
        error('Lipschitz constant L must be positive.');
    end
    if sigma < 0
        error('Strong convexity parameter sigma cannot be negative.');
    end
    % Prevent q >= 1 which happens if sigma >= L (or L is very small)
    if sigma >= L
       warning('AGD:SigmaGTE_L', 'Sigma >= L. AGD might not converge optimally. Setting q based on L/sigma.');
       % Avoid division by zero if sigma is also zero
       if sigma > 0
           q = L / sigma; % Use the inverse ratio if sigma >= L
       else
           q = 0; % Default to no strong convexity if sigma is zero or negative
           warning('AGD:SigmaZero', 'Sigma is non-positive. Setting q=0.');
       end
    else
       q = sigma/L; % Standard case
    end
    % Ensure q is within [0, 1) for stability of the momentum term
    q = max(0, min(q, 1 - 1e-9)); % Clamp q slightly below 1

    % Initialize array to store objective function values
    obj_values = zeros(max_iter, 1);
    f_smooth_values = zeros(max_iter, 1); % Internal storage for smoothed values
    dis_agd = zeros(max_iter, 1);
    convergence = false;
    % grad_norms = zeros(max_iter, 1);

    % Start timing
    tic;

    final_iter = 0; % Track actual iterations
    % --- AGM iterations ---
    for iter = 1:max_iter
        final_iter = iter;

        % Compute the objective function values (based on mu)
        obj = sum(exp(mu)) + sum(B .* max(log(v)-mu, [], 2)) - fval_solver;

        % Compute the smoothing function values (based on mu)
        max_log_v_mu = max(log(v) - repmat(mu, n, 1), [], 2);
        rescaled_log_v_mu = (log(v) - repmat(mu, n, 1) - max_log_v_mu) / delta;
        log_sum_exp_term = log(sum(exp(rescaled_log_v_mu), 2));
        log_sum_exp_term(isinf(log_sum_exp_term)) = -realmax; % Handle potential -Inf
        f_smooth = sum(exp(mu)) + delta * sum(B .* ((max_log_v_mu/delta) + log_sum_exp_term));

        % Document values for the current iteration 'iter' (based on state 'mu')
        obj_values(iter) = obj;
        f_smooth_values(iter) = f_smooth; % Store internal smoothed value
        dis_agd(iter) = norm(exp(mu) - p_opt_solver);

        % ! Stable gradient calculator (using y)
        temp1 = log(v) - repmat(y, n, 1); % y \in R^{1*m}
        max_temp1 = max(temp1,[],2); % Normalize for every row
        exp_temp1 = exp((temp1 - max_temp1)/delta);
        sum_exp_temp1 = sum(exp_temp1, 2);
        valid_rows = sum_exp_temp1 > 1e-100; % Use a small threshold for stability
        cal_temp1 = zeros(n, m);
        if any(valid_rows)
            cal_temp1(valid_rows, :) = exp_temp1(valid_rows, :) ./ sum_exp_temp1(valid_rows);
        end
        temp_2 = sum(B.*cal_temp1); % * n*1 by n*m
        grad_f = exp(y) - temp_2; % Gradient at y
        % ! Stable end

        % Store the gradient norm (optional)
        % grad_norms(iter) = norm(grad_f);

        % Update of mu and y
        step = 1 / L; % Step size based on L
        mu_new = P(y - step * grad_f);

        % Update y (Momentum step)
        momentum_factor = (1 - sqrt(q)) / (1 + sqrt(q));
        y_new = mu_new + momentum_factor * (mu_new - mu);

        % --- Stopping Criteria ---
        % 1. Primary Convergence Check (Original Objective Gap)
        if iter >= 2 && obj < epsilon
            convergence = true;
            break;
        end

        % 2. Adaptive Phase Termination Check (Smoothed Objective Plateau)
        if adaptive && iter >= 200 % Requires iter >= 200
            if abs(f_smooth_values(iter) - f_smooth_values(iter-1)) < 1e-3 && ...
               abs(f_smooth_values(iter-1) - f_smooth_values(iter-2)) < 1e-3 && ...
               abs(f_smooth_values(iter-2) - f_smooth_values(iter-3)) < 1e-3 && ...
               abs(f_smooth_values(iter-3) - f_smooth_values(iter-4)) < 1e-3
                fprintf('AGD Info: Adaptive phase termination triggered at iteration %d.\n', iter);
                break; % Exit loop (phase termination)
            end
        end
        % --- End Stopping Criteria ---

        % Update variables for next iteration
        mu = mu_new;
        y = y_new;

        % Check for NaN/Inf (optional safeguard)
        if any(isnan(mu)) || any(isinf(mu)) || any(isnan(y)) || any(isinf(y))
            warning('AGD:NaNInf', 'NaN or Inf detected in mu or y at iteration %d. Stopping.', iter);
            break; % Exit loop
        end

    end % End of iterations loop
    % --- End AGM iterations ---

    % End timing
    time = toc;

    % Extract the solution (use the last valid mu_new)
    solution = mu_new; % Or use 'mu' which corresponds to the last calculated objectives

    % Trim the arrays to the actual number of iterations performed
    obj_values = obj_values(1:final_iter);
    f_smooth_values_agd = f_smooth_values(1:final_iter); % Assign trimmed values to output
    dis_agd = dis_agd(1:final_iter);
    % grad_norms = grad_norms(1:iter);

    % Final convergence check
    if ~convergence && ~isempty(obj_values) && obj_values(end) < epsilon
        convergence = true;
    end
    if ~convergence && final_iter == max_iter
         fprintf('AGD Warning: Did not converge within %d iterations. Final objective gap: %e\n', max_iter, obj_values(end));
    end

    % --- Plotting ---
    % Plot the primary results if plot_flag is true
    if plot_flag
        figure; % Create figure for primary plots
        iterations = 1:final_iter;

        subplot(2, 1, 1);
        semilogy(iterations, abs(obj_values), 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('Original Objective Gap');
        title('AGD - Objective Value Convergence');
        grid on;
        hold on; yline(epsilon, 'r--', 'Tolerance'); hold off;

        subplot(2, 1, 2);
        plot(iterations, dis_agd, 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('Distance ||exp(mu) - p*||');
        title('AGD - Distance to Optimal Solution');
        grid on;
    end

    % Plot smoothed objective function value separately if plot_flag_smooth is true
    if plot_flag_smooth
        figure; % Create a new figure for the smooth plot
        iterations = 1:final_iter;
        plot(iterations, f_smooth_values_agd, 'LineWidth', 1.5); % Use the output variable
        xlabel('Iteration');
        ylabel('Smoothed Objective Value');
        title('AGD - Smoothed Objective Function Value');
        grid on;
    end
    % --- End Plotting ---

end % End of function
