function [solution, time, iter, obj_values, distance_md] = linear_primal_md_adaptive(v, B, x_0, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver, switch_step)
    % adaptive_linear_primal_md: Solves the optimization problem using Adaptive Mirror Descent.
    %
    % Inputs:
    % v: n x m matrix, given values
    % B: n x 1 vector, given constraints
    % x_0: n x m matrix, initial point
    % epsilon: stopping criterion
    % max_iter: maximum number of iterations
    % plot_flag: flag to plot results
    % p_opt_solver: 1 x m vector, optimal p values for comparison
    % fval_solver: optimal objective value for comparison
    % switch_step: number of iterations after which the step size decreases
    %
    % Outputs:
    % solution: final solution of the allocation problem
    % time: total computation time
    % iter: number of iterations performed
    % obj_values: array of objective function values over iterations
    % distance_md: array of distances to the optimal solution over iterations

    [n, m] = size(v);
    x = x_0; % Initialize x
    iter = 1;
    obj_values = []; % Array to store objective function values
    distance_md = []; % Array to store distances to the optimal solution

    tic;
    while iter < max_iter
        % Determine the current step size based on the iteration count
        if iter <= switch_step
            eta = 0.4; % Step size 4 for the first segment
        elseif iter <= 2 * switch_step
            eta = 0.45; % Step size 3 for the second segment
        elseif iter <= 3 * switch_step
            eta = 0.5; % Step size 2 for the third segment
        else
            eta = 0.5; % Step size 1 for the remaining iterations
        end

        % Compute gradient
        grad = -log(v) + log(sum(x, 1));

        % Compute current p_j = sum_i x_{i j}
        p_current = sum(x, 1);
        distance_current = norm(p_current - p_opt_solver, 2);
        distance_md = [distance_md, distance_current];

        % Compute objective value
        obj = sum(p_current) - sum(B .* log(min(p_current ./ v, [], 2))) - fval_solver;
        obj_values = [obj_values, obj];

        % Update step using exponential gradient descent
        x_temp = x .* exp(-eta * grad);
        sum_temp = repmat(sum(x_temp, 2) ./ B, 1, m);
        x_new = x_temp ./ sum_temp;

        % Stopping criteria
        if iter > 2 && obj < epsilon
            break;
        end

        % Update x
        x = x_new;
        iter = iter + 1;
    end
    time = toc;
    solution = sum(x, 1);

    % Plot the figure if plot_flag is 1
    if plot_flag
        figure;
        subplot(2, 1, 1);
        plot(1:length(obj_values), obj_values, '-o');
        xlabel('Iteration');
        ylabel('Function Value Gap');
        title('Adaptive MD - Function Value Convergence');

        subplot(2, 1, 2);
        plot(1:length(distance_md), distance_md, '-o');
        xlabel('Iteration');
        ylabel('Iteration Distance');
        title('Adaptive MD - Iteration Convergence');
    end
end