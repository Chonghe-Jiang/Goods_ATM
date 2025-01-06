function [solution, time, iter, obj_values, distance_md] = quasi_primal_md(v, B, x_0, eta, epsilon, max_iter, plot_flag, p_opt_solver, fval_solver)
    [n, m] = size(v);
    % mirror_descent: Solves the following optimization problem using Mirror Descent:
    %
    % max sum_{i, j} b_{i j} log(v_{i j}) - sum_j (sum_i b_{i j}) log(sum_i b_{i j})
    % s.t. sum_j b_{i j} = B_i, i in [n]
    %      b_{i j} >= 0
    %%% ! The above is the md for problem in linear case
    %%% ! Note that it is different from goods primal (we do some variable replacement)
    %%% ! Now let us share how to do variable recovery
    %%% Todo: check the correctness of the MD for linear case - B_i = 1 is correct
    %%% ! Yes, it is right
    %%% ! The below is the md for quasi linear case, see page 7 of first order paper
    % \min _{\bar{b}=(b, \delta)} \varphi(b)=-\sum_{i, j}\left(1+\log v_{i j}\right) 
    % b_{i j}+\sum_j p_j(b) \log p_j(b)
    %  \text { s.t. } \bar{b} \in \mathcal{B} 
    %%% Todo: what we need to do for quasi linear - check correctness similar to linear - change the gradient writing
    %%% ! function formulation
    %%% ! gradient
    %%% ! recover
    %%% ! initialization -> should refer to other file to check


    x_0_extended = x_0;
    x = x_0_extended;
    
    iter = 1;
    obj_values = []; % Array to store objective function values
    distance_md = [];
    tic;
    while iter < max_iter
        % Compute gradient
        % grad_{i, j} = 1 + log(v_{i, j}) - log(sum_i x_{i, j}) for j = 1:m
        % grad_{i, m+1} = 0 (since δ_i does not directly appear in the objective function)
        grad = zeros(n, m+1);
        grad(:, 1:m) = 1 + log(v) - log(sum(x(:, 1:m), 2));
        grad(:, m+1) = 0; % Gradient for δ_i is zero
        
        % Compute current p_j = sum_i x_{i j}
        p_current = sum(x(:, 1:m), 1);
        distance_current = norm(p_current - p_opt_solver, 2);
        distance_md = [distance_md, distance_current];
        
        % Compute objective function value
        obj = sum(p_current) - sum(B .* log(min(p_current ./ v, [], 2))) - fval_solver;
        obj_values = [obj_values, obj];
        
        % Update step using exponential gradient descent
        x_temp = x .* exp(-eta * grad);
        
        % Normalize to satisfy the constraint sum_j b_ij + δ_i = B_i
        sum_temp = repmat(sum(x_temp, 2) ./ B, 1, m+1);
        x_new = x_temp ./ sum_temp;
        
        % Stopping criteria by using the distance of the point p: gap(p)/stepsize
        if iter > 2 && obj < epsilon
            break;
        end

        % Update x
        x = x_new;
        iter = iter + 1;
    end
    time = toc;
    x_mirror = x;
    solution = sum(x_mirror(:, 1:m), 1);

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