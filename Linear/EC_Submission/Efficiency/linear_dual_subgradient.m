function [solution,obj_values, dis_sub, time, iter] = linear_dual_subgradient(v, B, p0, max_iter, step_size, epsilon, plot_flag, p_opt_solver, fval_solver)
    %%% Todo - check whether this stopping criteria is okay
    % Input:
    % v - parameter matrix v \in R^{n*m}
    % B - vector B \in R^{n*1}
    % p0 - initial variable p \in R^{1*m}
    % max_iter - maximum number of iterations
    % step_size - step size for the subgradient method
    % epsilon - tolerance for convergence
    % plot_flag - true to plot, false to not plot
    % p_opt_solver - optimal solution from the mosek solver

    % Output:
    % solution - solution of the dual subgradient method
    % time - time taken to solve the problem
    % obj_values - array of objective function values at each iteration
    % dis_sub - array of distances between each iteration point and p_opt_solver

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Initialize the variable
    p = p0;

    % Initialize array to store objective function values
    obj_values = zeros(max_iter, 1);

    % Initialize array to store distances
    dis_sub = zeros(max_iter, 1);

    % Start timing
    tic;

    % Subgradient method iterations
    for iter = 1:max_iter
        % Compute the objective function and subgradient
        obj = sum(p) - sum(B .* log(min(p ./ v, [], 2))) - fval_solver;
        subgrad = ones(1, m);
        for i = 1:n
            [min_val, min_idx] = min(p ./ v(i, :));
            selected_idx = min_idx(randi(length(min_idx)));
            subgrad(selected_idx) = subgrad(selected_idx) - B(i) / p(selected_idx);
        end

        % Store the objective function value
        obj_values(iter) = obj;

        % Calculate the distance to p_opt_solver
        dis_sub(iter) = norm(p - p_opt_solver);

        % Update the variable
        p_new = p - step_size * subgrad;

        % Project onto the feasible set (non-negative values)
        p_new = max(p_new, 0);

        % Check convergence - by distance of the p vector
        %%% Todo - check whether this stopping criteria is okay
        if iter >=2 && obj < epsilon
             break;
        end

        % Update the variable
        p = p_new;
    end

    % End timing
    time = toc;

    % Extract the solution
    solution = p;

    % Trim the obj_values array to the actual number of iterations
    obj_values = obj_values(1:iter);

    % Trim the dis_sub array to the actual number of iterations
    dis_sub = dis_sub(1:iter);

    % Plot the results if plot_flag is true
    if plot_flag
        figure;
        iterations = 1:iter;
        subplot(2, 1, 1);
        plot(iterations, abs(obj_values), 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Function Value Gap');
        title('Sub - Function Value Convergence');
        grid on;

        subplot(2, 1, 2);
        plot(iterations, dis_sub, 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Iteration Distance');
        title('Sub - Iteration Convergence');
        grid on;
    end
end