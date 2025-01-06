function [solution, time, obj_values, f_smooth_values, dis_agd] = linear_dual_adaptive_agd(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, p_opt_solver, fval_solver, phase, epsilon_phase_change)
    % Input:
    % v - parameter matrix v \in R^{n*m}
    % B - vector B \in R^{n*1}
    % mu_0 - initial point mu_0 \in R^{1*m}
    % max_iter - maximum number of iterations
    % L - Lipschitz constant
    % sigma - strong convexity parameter
    % epsilon - tolerance
    % mu_lower - lower bound for mu
    % mu_upper - upper bound for mu
    % delta - smoothing parameter
    % plot_flag - flag to plot results
    % p_opt_solver - optimal solution from another solver for distance calculation
    % fval_solver - optimal function value from another solver
    % phase - number of phases
    % epsilon_phase_change - tolerance for phase change
    % Output:
    % solution - solution of the optimization problem
    % time - time taken to solve the problem
    % obj_values - array of objective function values
    % f_smooth_values - array of smoothed objective function values
    % dis_agd - array of distances between each iteration point and p_opt_solver

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Define the projection operator
    P = @(mu) max(mu_lower, min(mu, mu_upper));

    % Initialize the variables and parameters
    mu = mu_0;
    y = mu_0;
    q = 0.001;

    % Initialize array to store objective function values
    obj_values = zeros(max_iter * phase, 1);
    f_smooth_values = zeros(max_iter * phase, 1);
    dis_agd = zeros(max_iter * phase, 1);

    % Start timing
    tic;

    % Total iterations counter
    total_iter = 0;

    % Array to store iterations per phase
    iter_per_phase = zeros(phase, 1);

    % Adaptive AGM phases
    for p = 1:phase
        % Update parameters for the current phase
        if p > 1
            delta = delta / 2;
            q = q / 2;
            L = L + sum(B) / delta;
        end

        % AGM iterations for the current phase
        for iter = 1:max_iter
            total_iter = total_iter + 1;

            % Compute the objective function values
            obj = sum(exp(mu)) + sum(B .* max(log(v)-mu, [], 2));
            f_smooth = sum(exp(mu)) + delta * sum(B .* log(sum(exp((log(v) - mu) / delta), 2)));
            obj_values(total_iter) = obj;
            f_smooth_values(total_iter) = f_smooth;

            % Calculate the distance to p_opt_solver
            dis_agd(total_iter) = norm(exp(mu) - p_opt_solver);

            % Check if the distance to p_opt_solver is less than epsilon
            if dis_agd(total_iter) < epsilon
                % End timing
                time = toc;

                % Extract the solution
                solution = mu;

                % Trim the obj_values array to the actual number of iterations
                obj_values = obj_values(1:total_iter);
                f_smooth_values = f_smooth_values(1:total_iter);
                dis_agd = dis_agd(1:total_iter);

                % Record iterations per phase
                iter_per_phase(p) = iter;

                % Print phase and iteration information
                fprintf('Total phases: %d\n', p);
                for i = 1:p
                    fprintf('Phase %d: %d iterations\n', i, iter_per_phase(i));
                end
                fprintf('Total iterations: %d\n', total_iter);

                % Plot the results if plot_flag is true
                if plot_flag
                    figure;
                    iterations = 1:total_iter;
                    subplot(2, 1, 1);
                    plot(iterations, abs(obj_values - fval_solver), 'LineWidth', 2);
                    xlabel('Iteration');
                    ylabel('Function Value Gap');
                    title('Function Value Gap');
                    grid on;

                    subplot(2, 1, 2);
                    plot(iterations, dis_agd, 'LineWidth', 2);
                    xlabel('Iteration');
                    ylabel('Distance to p_opt_solver');
                    title('Distance Convergence');
                    grid on;
                end

                return;
            end

            % Update mu in the current iteration
            % gradient of f_delta - need to normalize the variable mu
            temp1 = log(v) - repmat(mu, n, 1);
            max_temp1 = repmat(max(temp1, [], 2), 1, m);

            temp1 = (temp1 - max_temp1) / delta;  % enhance computation stability
            temp1 = exp(temp1); 
            temp2 = sum(temp1, 2) ./ (B);
            temp2 = repmat(temp2, 1, m);
            temp2 = temp1 ./ temp2;
            grad_f = exp(mu) - sum(temp2);
            % Update of mu and y
            mu_new = P(y - 1 / L * grad_f);
            % Update y
            y_new = mu_new + ((1 - sqrt(q)) / (1 + sqrt(q))) * (mu_new - mu);

            % Check changing phase phenomena
            if iter >= 2 && norm(grad_f) < epsilon_phase_change
                break;
            end

            % Update variables
            mu = mu_new;
            y = y_new;
        end

        % Record iterations per phase
        iter_per_phase(p) = iter;

        % Print grad_f norm at the end of the current phase
        fprintf('Phase %d ended with grad_f norm: %f\n', p, norm(grad_f));
    end

    % End timing
    time = toc;

    % Extract the solution
    solution = mu_new;

    % Trim the obj_values array to the actual number of iterations
    obj_values = obj_values(1:total_iter);
    f_smooth_values = f_smooth_values(1:total_iter);
    dis_agd = dis_agd(1:total_iter);

    % Print phase and iteration information
    fprintf('Total phases: %d\n', phase);
    for i = 1:phase
        fprintf('Phase %d: %d iterations\n', i, iter_per_phase(i));
    end
    fprintf('Total iterations: %d\n', total_iter);

    % Plot the results if plot_flag is true
    if plot_flag
        figure;
        iterations = 1:total_iter;
        subplot(2, 1, 1);
        plot(iterations, abs(obj_values - fval_solver), 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Function Value Gap');
        title('Function Value Gap');
        grid on;

        subplot(2, 1, 2);
        plot(iterations, dis_agd, 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Distance to p_opt_solver');
        title('Distance Convergence');
        grid on;
    end
end