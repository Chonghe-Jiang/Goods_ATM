function [solution, time, iter, obj_values, distance_da ] = linear_primal_da(v, B, x_0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver)
    % linear_primal_da: Solves the optimization problem using Dual Averaging.
    % Includes logging for the initial point (iteration 0).
    %
    % Inputs:
    % v: n x m matrix, given values (should be positive)
    % B: n x 1 vector, given constraints (row sums)
    % x_0: n x m matrix, initial point (should satisfy constraints or be projected)
    % eta: step size parameter (influences beta_t)
    % epsilon: stopping criterion threshold for objective gap
    % max_iter: maximum number of iterations
    % plot_flag: boolean, whether to plot convergence graphs
    % p_opt_solver: 1 x m vector, optimal p values (column sums) for distance comparison
    % fval_solver: scalar, optimal function value for gap calculation
    %
    % Outputs:
    % solution: 1 x m vector, final column sums (sum_i x_{i j})
    % time: scalar, execution time
    % iter: scalar, number of iterations performed (excluding initial point)
    % obj_values: array, history of objective function gap values (includes initial point)
    % distance_da: array, history of L2 distance to optimal p values (includes initial point)

    [n, m] = size(v);

    % Ensure v has positive values to avoid log(0) or log(-)
    if any(v(:) <= 0)
        error('Matrix v must contain positive values.');
    end

    % Initialize primal variable x
    x = x_0;

    % Check if initial x satisfies constraints (optional but recommended)
    initial_row_sums = sum(x, 2);
    if max(abs(initial_row_sums - B)) > 1e-6 % Tolerance for floating point errors
        warning('Initial point x_0 does not satisfy row sum constraints. Results may be inaccurate.');
        % Optional: Project x_0 onto the feasible set here if needed.
        % x = projection_simplex_rows(x_0, B); % Assuming such a function exists
    end

    % Initialize dual average accumulator
    z = zeros(n, m);

    % Initialize arrays to store objective function values (gap) and distances
    % Pre-allocate slightly larger if needed, or use dynamic arrays
    obj_values = [];
    distance_da = [];

    % --- Calculate and Log Values for Initial Point (Iteration 0) --- START REVISION
    p_initial = sum(x, 1); % p based on x_0
    p_initial(p_initial <= 0) = 1e-15; % Stability for log
    min_ratios_initial = min(p_initial ./ v, [], 2); % min along dim 2 (columns)
    min_ratios_initial(min_ratios_initial <= 0) = 1e-15; % Stability for log
    obj_initial = sum(p_initial) - sum(B .* log(min_ratios_initial)) - fval_solver;
    distance_initial = norm(p_initial - p_opt_solver, 2);

    obj_values = [obj_values, obj_initial]; % Log initial objective gap
    distance_da = [distance_da, distance_initial]; % Log initial distance
    % --- Calculate and Log Values for Initial Point (Iteration 0) --- END REVISION

    tic; % Start timer *after* initial point calculation

    iter = 0; % Initialize iteration counter (actual updates start from iter=1)
    while iter < max_iter % Loop for max_iter *updates*

        iter = iter + 1; % Increment iteration counter *before* calculations for this iter

        % --- Gradient Calculation ---
        % Compute current column sums p_j = sum_i x_{i j} (based on x from previous iter)
        p_current = sum(x, 1);
        p_current(p_current <= 0) = 1e-15; % Stability
        grad = log(p_current) - log(v); % Gradient based on x

        % --- Dual Accumulator Update ---
        z = z + grad;

        % --- Primal Variable Update (Solve DA Subproblem) ---
        beta_t = sqrt(iter) / eta; % Use current iteration number
        x_update_vals = exp(-z / beta_t);
        row_sums = sum(x_update_vals, 2);
        row_sums(row_sums <= 0) = 1e-15;
        x_new = x_update_vals .* (B ./ row_sums);

        % --- Logging for the *new* point x_new (Iteration 'iter') --- START REVISION
        p_new = sum(x_new, 1);
        p_new(p_new <= 0) = 1e-15; % Stability for log
        min_ratios_per_row = min(p_new ./ v, [], 2); % min along dim 2 (columns)
        min_ratios_per_row(min_ratios_per_row <= 0) = 1e-15; % Stability for log
        obj = sum(p_new) - sum(B .* log(min_ratios_per_row)) - fval_solver;
        distance_current = norm(p_new - p_opt_solver, 2);

        obj_values = [obj_values, obj]; % Append objective gap for iteration 'iter'
        distance_da = [distance_da, distance_current]; % Append distance for iteration 'iter'
         % --- Logging for the *new* point x_new (Iteration 'iter') --- END REVISION

        % --- Stopping criterion ---
        % Check based on the objective value calculated for the *current* iteration 'iter'
        if iter > 0 && abs(obj) < epsilon % Check absolute value of gap (iter > 0 check is implicit)
             break; % Exit loop if converged
        end

        % Update primal variable for the next iteration
        x = x_new;

    end % End of while loop

    time = toc; % Stop timer

    % Final solution (column sums from the last computed x)
    solution = sum(x, 1);

    % Plotting (if requested)
    if plot_flag
        figure; % Create a new figure for DA results
        plot_iterations = 0:iter; % Iterations from 0 (initial) to final iter

        subplot(2, 1, 1); % First subplot: Objective Function Gap
        plot(plot_iterations, obj_values, '-o'); % Plot includes initial point
        xlabel('Iteration');
        ylabel('Function Value Gap');
        title('DA - Function Value Convergence');
        grid on;

        subplot(2, 1, 2); % Second subplot: Iteration Distance
        plot(plot_iterations, distance_da, '-o'); % Plot includes initial point
        xlabel('Iteration');
        ylabel('L2 Distance to p_{opt}');
        title('DA - Iteration Convergence (Distance to p_{opt})');
        grid on;
    end

    % Return the actual number of *update* iterations performed
    % iter value already holds this count correctly

end