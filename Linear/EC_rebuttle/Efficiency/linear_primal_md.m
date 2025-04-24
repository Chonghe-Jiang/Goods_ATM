function [solution, time, iter, obj_values, distance_md ] = linear_primal_md(v, B, x_0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver)
    % mirror_descent: Solves the following optimization problem using Mirror Descent:
    %
    % max sum_{i, j} b_{i j} log(v_{i j}) - sum_j (sum_i b_{i j}) log(sum_i b_{i j})
    % s.t. sum_j b_{i j} = B_i, i in [n]
    %      b_{i j} >= 0
    
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
    % x_final: n x m matrix, final solution of the allocation problem

    [n, m] = size(v);
    
    % Project initial point x_0 to satisfy constraints
    % x = projection_simplex(x_0, B, n, m);
    %%% Todo - same initilization problem - solved
    x = x_0;
    iter = 1;
    obj_values = []; % Array to store objective function values
    distance_md = [];
    tic;
    while iter < max_iter
        % Compute gradient
        % grad_{i, j} = log(v_{i, j}) - log(sum_i x_{i, j}) - 1
        % grad = log(v) - log(sum(x, 1)) - 1; - wrong version
        % grad = 1 - log(v) + log(sum(x, 1));
        grad = - log(v) + log(sum(x, 1));
        % Todo- please check the gradient computation
        %%% Warning: This is a max problem so we should convert it into a min problem
        %%% Warning: Consider the gradient here - gradient descent of the negative objective function
        % Compute current p_j = sum_i x_{i j}
        p_current = sum(x, 1);
        distance_current = norm(p_current - p_opt_solver, 2);
        distance_md = [distance_md, distance_current];
        obj = sum(p_current) - sum(B .* log(min(p_current ./ v, [], 2))) - fval_solver; 
        obj_values = [obj_values, obj];
        % Update step using exponential gradient descent
        % x_new_{i, j} = x_{i, j} * exp(-eta * grad_{i, j})
        % Use the normalization from the Bolte paper (already theorem)
        x_temp = x .* exp(-eta * grad);
        sum_temp = repmat(sum(x_temp,2)./B,1,m); 
        x_new = x_temp./sum_temp;

        % Stopping criteria by using the distance of the point p: gap(p)/stepsize
        % Alternative1 - with the solver
        % Alternative2 - using the function value
        %%% Todo - unifying the stopping creteria
        if iter > 2 && obj < epsilon
            break;
        end

        % Update x
        x = x_new;
        iter = iter + 1;
    end
    time = toc;
    x_mirror = x;
    solution = sum(x_mirror, 1);


    % Plot the figure if plot_figure is 1
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

