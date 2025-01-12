function [solution, time, iter, obj_values, dis_agd, convergence] = quasi_dual_agd(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive)
    %%% ! we change the objective function calculation and the gradient calculation

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Define the projection operator
    P = @(mu) max(mu_lower, min(mu, mu_upper));

    % Initialize the variables and parameters
    %%% Todo - Initialization problem
    mu = mu_0;
    y = mu_0;
    % q = 0.1; 
    %%% Todo - choose the best parameter here
    q = sigma/L; % L is only multiplied here

    % Initialize array to store objective function values
    obj_values = zeros(max_iter, 1);
    f_smooth_values = zeros(max_iter, 1);
    dis_agd = zeros(max_iter, 1);
    convergence = false;
    % grad_norms = zeros(max_iter, 1);

    % Start timing
    tic;

    % AGM iterations
    for iter = 1:max_iter
        % Compute the objective function values
        obj = sum(exp(mu)) + sum(B .* max(log(v)-mu, [], 2)) - fval_solver;

        % Compute the smoothing function values
        % ! Unstable: f_smooth = sum(exp(mu)) + delta * sum(B .* log(sum(exp((log(v) - repmat(mu,n,1)) / delta), 2))); 
        % ! Solved - New stable one - minus the maximizer
        % * This resolve the function problem
        max_log_v_mu = max(log(v) - repmat(mu, n, 1), [], 2);
        % Rescale the values by subtracting the minimum
        rescaled_log_v_mu = (log(v) - repmat(mu, n, 1) - max_log_v_mu) / delta;

        % Compute the log-sum-exp term using the rescaled values
        log_sum_exp_term = log(1 + sum(exp(rescaled_log_v_mu), 2));

        % Compute the final expression
        f_smooth = sum(exp(mu)) + delta * sum(B .* ((max_log_v_mu/delta) + log_sum_exp_term));
        % ! New stable end


        % Document some values
        obj_values(iter) = obj;
        f_smooth_values(iter) = f_smooth; 
        dis_agd(iter) = norm(exp(mu) - p_opt_solver); 
        % ! Stable gradient calculator
        temp1 = log(v) - repmat(y, n, 1); % y \in R^{1*m}
        max_temp1 = max(temp1,[],2); % Normalize for every row
        exp_temp1 = exp((temp1 - max_temp1)/delta);
        cal_temp1 = exp_temp1./ (1 + sum(exp_temp1,2));
        temp_2 = sum(B.*cal_temp1); % * n*1 by n*m
        grad_f = exp(y) - temp_2;
        % ! Stable end
        
        % Store the gradient norm
        % grad_norms(iter) = norm(grad_f);

        % Update of mu and y
        mu_new = P(y - (1 / (2*L)) * grad_f);
        
        % Update y
        y_new = mu_new + ((1 - sqrt(q)) / (1 + sqrt(q))) * (mu_new - mu);
        
        % Check convergence
        %%% Todo: how the stepsize is used in the stopping criteria in the three algorithms -whether we can use this information or not
        if iter >= 2 && obj < epsilon
            convergence = true;
            break;
        end
        %%% Todo: give a rigorous definition of the phase changing phenomena - several selections: gradient, relative gradient, gap  
        if adaptive && iter>=20 && abs(f_smooth_values(iter) - f_smooth_values(iter-1)) < 1e-3 &&  abs(f_smooth_values(iter-1) - f_smooth_values(iter-2)) < 1e-3 && abs(f_smooth_values(iter-2) - f_smooth_values(iter-3)) < 1e-3 &&  abs(f_smooth_values(iter-3) - f_smooth_values(iter-4)) < 1e-3
            break;
        end
        
        % Update variables
        mu = mu_new;
        y = y_new;
    end

    % End timing
    time = toc;

    % Extract the solution
    solution = mu_new;

    % Trim the obj_values array to the actual number of iterations
    obj_values = obj_values(1:iter);
    % f_smooth_values = f_smooth_values(1:iter);
    dis_agd = dis_agd(1:iter);
    % grad_norms = grad_norms(1:iter);

    % Plot the results if plot_flag is true
    if plot_flag
        figure;
        iterations = 1:iter;
        
        subplot(2, 1, 1);
        plot(iterations, abs(obj_values), 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Function Value Gap');
        title('SAG - Function Value Convergence');
        grid on;
        
        subplot(2, 1, 2);       
        plot(iterations, dis_agd, 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Iteration Distance');
        title('SAG - Iteration Convergence');
        grid on;
        %%% Todo - Now we do not have an output for the smoothing function value
        % subplot(3, 1, 3);
        % plot(iterations, f_smooth_values, 'LineWidth', 2);
        % xlabel('Iteration');
        % ylabel('Smoothed Function Value');
        % title('Smoothed Function Value Convergence');
        % grid on;
    end
    if plot_flag_smooth
        figure;
        iterations = 1:iter;
        plot(iterations,f_smooth_values(1:iter),'LineWidth',2);
        xlabel('Iteration');
        ylabel('Smoothing Function Value');
        title('Smoothing Function Value');
    end
end