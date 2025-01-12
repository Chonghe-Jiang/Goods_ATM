function [solution, time, iter, obj_values, dis_agd] = linear_dual_rounding(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, p_opt_solver, fval_solver)
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
    % p_opt_solver - optimal solution from another solver for distance calculation
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
    q = sigma/L;

    % Initialize array to store objective function values
    obj_values = zeros(max_iter, 1);
    dis_agd = zeros(max_iter, 1);
    % Start timing
    tic;

    % AGM iterations
    for iter = 1:max_iter
        % Compute the objective function values
        obj = sum(exp(mu)) + sum(B .* max(log(v)-mu, [], 2));
        % Store the objective function value
        obj_values(iter) = obj;
        % Calculate the distance to p_opt_solver
        dis_agd(iter) = norm(exp(mu) - p_opt_solver);

        % Update mu in the current iteration
        % gradient of f_delta - need to normalize the variable mu
        temp1 = log(v) - repmat(y, n, 1); % y \in R^{1*m}
        max_temp1 = max(max(temp1));

        temp1 = (temp1 - max_temp1) / delta;  % enhance computation stability
        temp1 = exp(temp1); % upper number
        temp2 = sum(temp1, 2) ./ (B);
        temp2 = repmat(temp2, 1, m);
        temp2 = temp1 ./ temp2;
        grad_f = exp(y) - sum(temp2);

        % Update of mu and y
        mu_new = P(y - (1 / (2*L)) * grad_f);
        
        % Update y
        y_new = mu_new + ((1 - sqrt(q)) / (1 + sqrt(q))) * (mu_new - mu);
        
        % Check convergence
        if iter >= 2 && abs(obj_values(iter)-fval_solver) < epsilon
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
end