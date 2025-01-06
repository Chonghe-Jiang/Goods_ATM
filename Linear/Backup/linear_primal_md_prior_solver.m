function [p_mirror, obj_mirror, obj_values] = linear_primal_md(v, B, x_0, eta, epsilon, max_iter, p_opt_solver, plot_flag)
    % mirror_descent: Solves the following optimization problem using Mirror Descent:
    %
    % max sum_{i, j} b_{i j} log(v_{i j}) - sum_j (sum_i b_{i j}) log(sum_i b_{i j})
    % s.t. sum_j b_{i j} = B_i, i in [n]
    %      b_{i j} >= 0
    %
    % Inputs:
    % v: n x m matrix, given values
    % B: n x 1 vector, given constraints
    % x_0: n x m matrix, initial point (will be projected to satisfy constraints)
    % eta: step size sequence
    % epsilon: stopping criterion
    % max_iter: maximum number of iterations
    % p_opt_solver: 1 x m vector, optimal p values for comparison
    % plot_flag
    %
    % Outputs:
    % x_final: n x m matrix, final solution of the allocation problem

    [n, m] = size(v);
    
    % Project initial point x_0 to satisfy constraints
    % x = projection_simplex(x_0, B, n, m);
    x = x_0;
    iter = 0;
    % distances = [];
    obj_values = []; % Array to store objective function values

    while iter < 10000
        % Compute gradient
        % grad_{i, j} = log(v_{i, j}) - log(sum_i x_{i, j}) - 1
        % grad = log(v) - log(sum(x, 1)) - 1; - wrong version
        grad = 1 - log(v) + log(sum(x, 1));
        % Warning: This is a max problem so we should convert it into a min problem
        p_current = sum(x, 1);
        obj = sum(p_current) - sum(B .* log(min(p_current ./ v, [], 2))); 
        obj_values = [obj_values, obj];
        % Update step using exponential gradient descent
        % x_new_{i, j} = x_{i, j} * exp(-eta * grad_{i, j})
        % Use the normalization from the Bolte paper (already theorem)
        x_temp = x .* exp(-eta * grad);
        sum_temp = repmat(sum(x_temp,2)./B,1,m);
        x_new = x_temp./sum_temp;

        % Compute current p_j = sum_i x_{i j}
        p_current = sum(x_new, 1);
        

        % % Compute distance to p_opt_solver
        % distance = norm(p_current - p_opt_solver, 2);
        % distances = [distances, distance];

        % % Compute objective function value
        % obj_value = sum(sum(x_new .* log(v))) - sum(p_current .* log(p_current));
        % objective_values = [objective_values, obj_value];

        % Check stopping criterion
        % If the change in x is smaller than epsilon, stop the iteration
        % if norm(x_new - x, 'fro') < epsilon
        %     break;
        % end

        % Update x
        x = x_new;
        iter = iter + 1;
    end

    x_mirror = x;
    p_mirror = sum(x_mirror, 1);
    obj_mirror = sum(p_mirror) - sum(B .* log(min(p_mirror ./ v, [], 2)));


    % Plot the figure if plot_figure is 1
    % if plot_flag
    %     figure;
    %     subplot(2, 1, 1);
    %     plot(1:length(distances), distances, '-o');
    %     xlabel('Iteration');
    %     ylabel('Distance to optimal');
    %     title('Distance Convergence');

    %     subplot(2, 1, 2);
    %     plot(1:length(objective_values), objective_values, '-o');
    %     xlabel('Iteration');
    %     ylabel('Objective Function Value');
    %     title('Function Convergence');
    % end
end

