function [solution, time, iter, obj_values, distance_md] = quasi_primal_md_adaptive(v, B, x_0, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver, switch_step)
    % mirror_descent: Solves the following optimization problem using Mirror Descent:
    %%% ! We are solving quasi-linear problem
    % max sum_{i, j} b_{i j} (1 + log(v_{i j})) - sum_j (sum_i b_{i j}) log(sum_i b_{i j})
    % s.t. sum_j b_{i j} + delta_i = B_i, i in [n]
    %      b_{i j} >= 0
    %      delta_i >= 0

    % Inputs:
    % v: n x m matrix, given values
    % B: n x 1 vector, given constraints
    % x_0: n x m matrix, initial point (will be projected to satisfy constraints)
    % eta: step size sequence
    % epsilon: stopping criterion
    % max_iter: maximum number of iterations
    % p_opt_solver: 1 x m vector, optimal p values for comparison
    % plot_flag

    % Outputs:
    % solution: 1 x m vector, final solution of the allocation problem
    % time: total computation time
    % iter: number of iterations
    % obj_values: array of objective function values
    % distance_md: array of distances to the optimal solution

    [n, m] = size(v);
    %%% ! Todo: fix the initialization issue -> still need to check
    % Initialize variables
    %%% ! Todo: this is also a change
    %%% ! X is the matrix with big size
    X = [x_0, B - sum(x_0, 2)]; % Augmented matrix: [x, delta]
    iter = 1;
    obj_values = []; % Array to store objective function values
    distance_md = []; % Array to store distances to the optimal solution
    tic;

    while iter < max_iter
        if iter <= switch_step
            eta = 0.3; % Step size 4 for the first segment
        elseif iter <= 2 * switch_step
            eta = 0.35; % Step size 3 for the second segment
        elseif iter <= 3 * switch_step
            eta = 0.4; % Step size 2 for the third segment
        else
            eta = 0.5; % Step size 1 for the remaining iterations
        end
        % Compute gradient
        % grad_{i, j} = 1 + log(v_{i, j}) - log(sum_i x_{i, j}) - 1
        % Simplify: grad_{i, j} = log(v_{i, j}) - log(sum_i x_{i, j})
        grad_x = -log(v) + log(sum(X(:, 1:m), 1)); % ! Change for gradient

        % Gradient for delta_i is 0 because delta_i does not appear in the objective
        grad_delta = zeros(n, 1);

        % Combine gradients into a single matrix
        grad = [grad_x, grad_delta];

        % Compute current p_j = sum_i x_{i j}
        p_current = sum(X(:, 1:m), 1);
        distance_current = norm(p_current - p_opt_solver, 2);
        distance_md = [distance_md, distance_current];

        % Compute objective function value
        %%% ! Here we need to change the expression of the original function
        obj = sum(p_current) - sum(B .* log(min([ones(n, 1), p_current ./ v], [], 2))) - fval_solver; 
        obj_values = [obj_values, obj];
       
        % Update step using exponential gradient descent
        X_temp = X .* exp(-eta * grad);
        sum_temp = repmat(sum(X_temp,2)./B,1,m+1); % ! Here we need to pay attention to the size
        % Normalize each row of X_temp to satisfy sum_j X_{i j} = B_i
        X_new = X_temp ./ sum_temp;

        % Stopping criteria
        if iter > 2 && obj < epsilon
            break;
        end

        % Update X
        X = X_new;
        iter = iter + 1;
    end

    time = toc;
    solution = sum(X(:, 1:m), 1);

    % Plot the figure if plot_flag is 1
    if plot_flag
        figure;
        subplot(2, 1, 1);
        plot(1:length(obj_values), obj_values, '-o');
        xlabel('Iteration');
        ylabel('Function Value Gap');
        title('MD - Function Value Convergence');

        subplot(2, 1, 2);
        plot(1:length(distance_md), distance_md, '-o');
        xlabel('Iteration');
        ylabel('Iteration Distance');
        title('MD - Iteration Convergence');
    end
end