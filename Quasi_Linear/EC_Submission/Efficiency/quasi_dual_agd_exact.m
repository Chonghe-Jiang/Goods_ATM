function [solution, time, iter, obj_values, dis_agd, convergence] = quasi_dual_agd_exact(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive)
    % Input:
    % v - parameter matrix v \in R^{n*m}
    % B - vector B \in R^{n*1}
    % mu_0 - initial point mu_0 \in R^{1*m}
    % L - Lipschitz constant
    % sigma - strong convexity parameter
    % epsilon - tolerance
    % mu_lower - lower bound for mu
    % mu_upper - upper bound for mu
    % delta - smoothing parameter
    % plot_flag - flag to plot results
    % p_opt_solver - optimal solution from the mosek solver
    % fval_solver - optimal function value from the mosek solver

    % Output:
    % solution - solution of the optimization problem
    % time - time taken to solve the problem
    % obj_values - array of objective function values
    % f_smooth_values - array of smoothed objective function values
    % dis_agd - array of distances between each iteration point and p_opt_solver
    % grad_norms - array of gradient norms at each iteration

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Define the projection operator
    P = @(mu) max(mu_lower, min(mu, mu_upper));

    % Initialize the variables and parameters
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
        %%% ! New gradient calculation - quasi-linear version
        %%% ! Change in function value
        max_log_v_mu = max([zeros(n, 1), log(v) - repmat(mu, n, 1)], [], 2); % Include 0 in the max
        obj = sum(exp(mu)) + sum(B .* max_log_v_mu) - fval_solver;

        % Compute the smoothing function values
        % Include 0 in the max and rescale
        temp1 = [zeros(n, 1), log(v) - repmat(mu, n, 1)]; % Include 0 in the max
        max_temp1 = max(temp1, [], 2); % Normalize for every row
        rescaled_temp1 = (temp1 - max_temp1) / delta; % Rescale
        exp_temp1 = exp(rescaled_temp1); % Exponentiate
        log_sum_exp_term = log(sum(exp_temp1, 2)); % Log-sum-exp

        % Compute the final smoothing function
        %%% ! Changes in calculatting the smooth function value
        f_smooth = sum(exp(mu)) + delta * sum(B .* ((max_temp1 / delta) + log_sum_exp_term));

        % Document some values
        obj_values(iter) = obj;
        f_smooth_values(iter) = f_smooth;
        %%% Todo: if use this distance, we should unify log or not
        dis_agd(iter) = norm(exp(mu) - p_opt_solver);

        % Compute the gradient
        % Include 0 in the softmax-like calculation
        % Stable gradient calculation
        % Include 0 in the stabilization
        %%% ! Calculating the gradient of y
        temp1 = [zeros(n, 1), log(v) - repmat(y, n, 1)]; % Include 0 in the max
        max_temp1 = max(temp1, [], 2); % Normalize for every row
        exp_temp1 = exp((temp1 - max_temp1) / delta); % Stabilized exponentials
        cal_temp1 = exp_temp1 ./ sum(exp_temp1, 2); % Softmax-like term, excluding the 0 term
        temp_2 = sum(B .* cal_temp1(:, 2:end)); % Exclude the 0 term for gradient
        grad_f = exp(y) - temp_2; % Gradient

        % Update of mu and y
        %%% Todo: Be careful here, whether we use long step or not
        mu_new = P(y - (4 / (L)) * grad_f); %%% Todo: previously 1/2L not 2/L - for synthetic data
        
        % Update y
        y_new = mu_new + ((1 - sqrt(q)) / (1 + sqrt(q))) * (mu_new - mu);
        
        %%% Todo: give a rigorous definition of the phase changing phenomena 
        if adaptive && iter>=500 && abs(f_smooth_values(iter) - f_smooth_values(iter-1)) < 1e-3 &&  abs(f_smooth_values(iter-1) - f_smooth_values(iter-2)) < 1e-3 && abs(f_smooth_values(iter-2) - f_smooth_values(iter-3)) < 1e-3 &&  abs(f_smooth_values(iter-3) - f_smooth_values(iter-4)) < 1e-3
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
    end
    %%% Current Version - No Smoothing Output
    if plot_flag_smooth
        figure;
        iterations = 1:iter;
        plot(iterations,f_smooth_values(1:iter),'LineWidth',2);
        xlabel('Iteration');
        ylabel('Smoothing Function Value');
        title('Smoothing Function Value');
    end
end