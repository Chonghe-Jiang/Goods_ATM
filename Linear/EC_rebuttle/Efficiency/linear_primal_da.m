function [solution, time, iter, obj_values, distance_da ] = linear_primal_da(v, B, x_0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver)
    % linear_primal_da: Solves the following optimization problem using Dual Averaging:
    %
    % min sum_j (sum_i x_{i j}) log(sum_i x_{i j}) - sum_{i, j} x_{i j} log(v_{i j})
    % s.t. sum_j x_{i j} = B_i, i in [n]
    %      x_{i j} >= 0
    %
    % This is equivalent to maximizing:
    % max sum_{i, j} x_{i j} log(v_{i j}) - sum_j (sum_i x_{i j}) log(sum_i x_{i j})
    
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
    
    % Outputs:
    % solution: 1 x m vector, final column sums (sum_i x_{i j})
    % time: scalar, execution time
    % iter: scalar, number of iterations performed
    % obj_values: array, history of objective function gap values
    % distance_da: array, history of L2 distance to optimal p values

    [n, m] = size(v);
    
    % Ensure v has positive values to avoid log(0) or log(-)
    if any(v(:) <= 0)
        error('Matrix v must contain positive values.');
    end
    
    % Initialize primal variable x
    % Option 1: Use the provided x_0 (assuming it's feasible)
    x = x_0; 
    % Option 2: Use a default feasible initialization (e.g., uniform)
    % x = repmat(B ./ m, 1, m); 
    
    % Check if initial x satisfies constraints (optional but recommended)
    initial_row_sums = sum(x, 2);
    if max(abs(initial_row_sums - B)) > 1e-6 % Tolerance for floating point errors
        warning('Initial point x_0 does not satisfy row sum constraints. Results may be inaccurate.');
        % Optional: Project x_0 onto the feasible set here if needed.
        % x = projection_simplex_rows(x_0, B); % Assuming such a function exists
    end

    % Initialize dual average accumulator
    z = zeros(n, m); 
    
    iter = 1;
    obj_values = []; % Array to store objective function values (gap)
    distance_da = []; % Array to store distance to optimal p

    tic; % Start timer
    
    while iter <= max_iter % Use <= for consistency with loop counter
        
        % --- Gradient Calculation ---
        % Compute current column sums p_j = sum_i x_{i j}
        p_current = sum(x, 1);
        
        % Avoid log(0) if p_current can be zero (unlikely if B > 0, v > 0)
        p_current(p_current <= 0) = 1e-15; % Add small epsilon for numerical stability
        
        % Gradient of g(x) = sum_j p_j log(p_j) - sum_{i,j} x_{ij} log(v_{ij})
        % grad_g_{i, j} = log(p_j) + 1 - log(v_{ij})
        % We use the gradient as computed in the MD code for consistency: log(p_j) - log(v_{ij})
        % This corresponds to minimizing the negative objective used in MD comments
        grad = log(p_current) - log(v); % Uses MATLAB's broadcasting

        % --- Dual Accumulator Update ---
        z = z + grad;
        
        % --- Primal Variable Update (Solve DA Subproblem) ---
        % Choose beta_t schedule. Example: beta_t = t / eta
        beta_t = sqrt(iter) / eta; 
        
        % Compute x_{t+1} = argmin_{x in X} { <z_{t+1}, x> + beta_{t+1} * psi(x) }
        % Using psi(x) = sum_{i,j} x_{ij} log(x_{ij}) (negative entropy)
        % Solution form: x_{t+1, ij} \propto exp(-z_{t+1, ij} / beta_{t+1})
        
        x_update_vals = exp(-z / beta_t);
        
        % Normalize rows to satisfy constraints sum_j x_{ij} = B_i
        row_sums = sum(x_update_vals, 2);
        
        % Avoid division by zero if a row_sum is zero (shouldn't happen if B_i > 0)
        row_sums(row_sums <= 0) = 1e-15; 
        
        % Apply normalization factor B_i / row_sum_i to each element in the row
        x_new = x_update_vals .* (B ./ row_sums); % Element-wise multiplication with broadcasting

        % --- Convergence Check & Logging ---
        % Calculate current objective gap (using the same formula as MD code)
        p_new = sum(x_new, 1);
        p_new(p_new <= 0) = 1e-15; % Stability for log
        % Calculate term based on dual objective or KKT conditions? 
        % The obj calculation seems specific, replicating it from MD code:
        % It looks potentially related to min(p_current ./ v, [], 2) which finds the min ratio per row.
        % Let's recalculate p_current based on the *previous* iterate x for the obj calculation,
        % or use p_new based on x_new. Let's use p_new for consistency with iteration t+1.
        
        % Calculate min ratio for each row i: min_j (p_new_j / v_ij)
        min_ratios_per_row = min(p_new ./ v, [], 2); % min along dim 2 (columns)
        min_ratios_per_row(min_ratios_per_row <= 0) = 1e-15; % Stability for log
        
        % Objective gap calculation (replicating MD code's logic)
        obj = sum(p_new) - sum(B .* log(min_ratios_per_row)) - fval_solver; 
        obj_values = [obj_values, obj]; % Append objective gap
        
        % Calculate distance to optimal p
        distance_current = norm(p_new - p_opt_solver, 2);
        distance_da = [distance_da, distance_current]; % Append distance

        % Stopping criterion (using the same logic as MD code)
        if iter > 1 && abs(obj) < epsilon % Check absolute value of gap
             break; % Exit loop if converged
        end

        % Update primal variable for the next iteration
        x = x_new;
        
        % Increment iteration counter
        iter = iter + 1;
    end
    
    time = toc; % Stop timer
    
    % Final solution (column sums)
    solution = sum(x, 1);

    % Plotting (if requested)
    if plot_flag
        figure; % Create a new figure for DA results
        
        subplot(2, 1, 1); % First subplot: Objective Function Gap
        plot(1:length(obj_values), obj_values, '-o');
        xlabel('Iteration');
        ylabel('Function Value Gap');
        title('DA - Function Value Convergence');
        grid on;

        subplot(2, 1, 2); % Second subplot: Iteration Distance
        plot(1:length(distance_da), distance_da, '-o');
        xlabel('Iteration');
        ylabel('L2 Distance to p_{opt}');
        title('DA - Iteration Convergence (Distance to p_{opt})');
        grid on;
    end
    
    % Ensure iter reflects the number of iterations completed
    if iter > max_iter 
        iter = max_iter; 
    end

end
