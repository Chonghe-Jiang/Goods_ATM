function [solution, obj_values, dis_sub, time, iter] = quasi_dual_subgradient(v, B, p0, max_iter, step_size, epsilon, plot_flag, p_opt_solver, fval_solver)
    %%% ! We change the nonsmooth subgradient calculation -> just an approximation one

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
        obj = sum(p) - sum(B .* log(max(p ./ v, 0, [], 2))) - fval_solver;
        subgrad = ones(1, m);
        for i = 1:n
            [max_val, max_idx] = max(p ./ v(i, :));
            if max_val >= 0
                selected_idx = max_idx(randi(length(max_idx)));
                subgrad(selected_idx) = subgrad(selected_idx) - B(i) / p(selected_idx);
            end
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
        if iter >= 2 && obj < epsilon
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