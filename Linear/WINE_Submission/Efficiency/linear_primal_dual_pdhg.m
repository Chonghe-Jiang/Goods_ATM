function [solution_x, solution_t, solution_p, solution_y, obj_values, dis_p, time, iter_total] = linear_primal_dual_pdhg(v, B, x0, t0, p0, y0, max_outer_iter, K_inner_iter, sigma_step, tau_step, epsilon, plot_flag, p_opt_solver, fval_solver)
    % linear_primal_dual_pdhg - Implements Primal-Dual Hybrid Gradient (PDHG) algorithm for Fisher Market Equilibrium.
    % Based on Algorithm 1 and equations (5), (6), (7) from "PDHCG: A Scalable First-Order Method..." (arXiv:2506.06258v1), Section 2.1.
    % This version explicitly uses auxiliary variable t and dual variable y.

    % Input:
    % v              - Valuation matrix u_ij (n x m)
    % B              - Budget vector w_i (n x 1)
    % x0             - Initial allocation matrix x (n x m)
    % t0             - Initial auxiliary variable t (n x 1)
    % p0             - Initial price vector p (1 x m)
    % y0             - Initial dual variable y (n x 1)
    % max_outer_iter - Maximum number of outer loop iterations (restarts)
    % K_inner_iter   - Number of inner iterations per outer loop
    % sigma_step     - Dual step size (sigma_k)
    % tau_step       - Primal step size (tau_k)
    % epsilon        - Tolerance for convergence (based on dual objective gap)
    % plot_flag      - True to plot, false to not plot
    % p_opt_solver   - Optimal price solution from external solver
    % fval_solver    - Optimal primal objective value from external solver

    % Output:
    % solution_x     - Final allocation matrix x
    % solution_t     - Final auxiliary variable t
    % solution_p     - Final price vector p
    % solution_y     - Final dual variable y
    % obj_values     - Array of dual objective function values (gap to optimal)
    % dis_p          - Array of distances between current p and p_opt_solver
    % time           - Time taken to solve the problem
    % iter_total     - Total number of inner iterations performed

    [n, m] = size(v);

    % Initialize variables for outer loop (Algorithm 1)
    x_n0 = x0; % x^(n,0)
    t_n0 = t0; % t^(n,0)
    p_n0 = p0; % p^(n,0)
    y_n0 = y0; % y^(n,0)
    
    % Ensure initial values are positive or within reasonable bounds
    % x_n0 = max(x_n0, 1e-12); % x >= 0 constraint
    % t_n0 = max(t_n0, 1e-12); % t > 0 constraint
    % p_n0 = max(p_n0, 1e-12); % p >= 0 constraint
    % y_n0 = y_n0; % y can be anything for now, check paper on y bounds


    % Store results (objective gap and price distance)
    % Allocate space for one extra initial point + max_outer_iter * K_inner_iter
    obj_values = zeros(max_outer_iter * K_inner_iter + 1, 1);
    dis_p = zeros(max_outer_iter * K_inner_iter + 1, 1);
    iter_idx = 1; % Start from 1 for the initial point

    tic; % Start overall timer

    % --- Record Initial Objective Function Value and Distance ---
    % Calculate objective value using initial p0 (p_n0)
    p_over_v_initial = p_n0 ./ v; % Corrected: Use p_n0 for initial calculation
    p_over_v_initial(p_over_v_initial <= 0) = 1e-12; % Protect against non-positive values
    
    beta_i_inv_initial = min(p_over_v_initial, [], 2);
    beta_i_inv_initial(beta_i_inv_initial <= 0) = 1e-12; % Ensure positivity for log

    initial_dual_obj = sum(p_n0) - sum(B .* log(beta_i_inv_initial)); % Corrected: Use p_n0
    obj_values(iter_idx) = abs(initial_dual_obj - fval_solver); % Gap to optimal dual value

    % Calculate distance to optimal prices using initial p0 (p_n0)
    dis_p(iter_idx) = norm(p_n0 - p_opt_solver); % Corrected: Use p_n0
    % -----------------------------------------------------------

    for outer_n = 0:max_outer_iter-1 % Outer loop: Restart mechanism
        % Initialize inner loop (Algorithm 1, Line 2)
        % These variables now hold the state *before* the first inner iteration update
        x_bar_nk = x_n0; % x_bar^(n,0)
        t_bar_nk = t_n0; % t_bar^(n,0)
        p_bar_nk = p_n0; % p_bar^(n,0)
        y_bar_nk = y_n0; % y_bar^(n,0)

        % Keep track of previous iterates for (2*var_k - var_prev) terms
        % For k=0 (first actual update), use var_n0 as prev values.
        x_prev = x_n0;
        t_prev = t_n0;
        p_prev = p_n0;
        y_prev = y_n0;

        % Set current iterates for the first inner loop update (k=0)
        x_k = x_n0;
        t_k = t_n0;
        p_k = p_n0;
        y_k = y_n0;

        for k = 0:K_inner_iter-1 % Inner loop. k=0 corresponds to the first actual update.
            iter_idx = iter_idx + 1; % Increment for the current iteration's results

            % Extrapolated primal variables for dual update (2*var_k - var_prev)
            x_extrapolated = 2 * x_k - x_prev;
            t_extrapolated = 2 * t_k - t_prev;

            % -------------------- Dual Variables Update (p, y) --------------------
            % (Algorithm 1, Line 3) and Eq. (5)
            
            % Update p (Equation 5)
            p_new = p_k + sigma_step * (sum(x_extrapolated, 1) - 1); % sum(x,1) sums columns
            % No explicit projection needed for p, but can add max(0) for safety if p starts negative
            % p_new = max(p_new, 1e-12); % Ensures non-negativity for later calculations

            % Update y (Equation 5)
            % sum_{j \in M} u_ij * x_ij is sum(v(i,:) .* x_extrapolated(i,:))
            y_new = y_k + sigma_step * (t_extrapolated - sum(v .* x_extrapolated, 2));
            % y can be negative, no explicit projection mentioned for y.

            % -------------------- Primal Variables Update (x, t) --------------------
            % (Algorithm 1, Line 4) and Eq. (6), (7)

            % Update x (Equation 7)
            % This is a simple projection onto R+ for each element
            x_new_unprojected = x_k - tau_step * (p_new - v .* y_new); % Element-wise operation
            x_new = max(x_new_unprojected, 0); % Proj_R+

            % Update t (Equation 6)
            % Ensure the term under sqrt is non-negative and B(i) is positive.
            sqrt_term_t = (tau_step * y_new - t_k).^2 + 4 * tau_step * B;
            t_new = (t_k - tau_step * y_new + sqrt(sqrt_term_t)) / 2;
            t_new = max(t_new, 1e-12); % Ensure t > 0 as per problem (2)

            % Store current and previous iterates for next iteration's extrapolation
            x_prev = x_k;
            t_prev = t_k;
            p_prev = p_k;
            y_prev = y_k;
            x_k = x_new;
            t_k = t_new;
            p_k = p_new;
            y_k = y_new;

            % -------------------- Update Averaged Solutions --------------------
            % (Algorithm 1, Line 5)
            % The averaging should now reflect k, which starts from 0 for the first update
            x_bar_nk = (k / (k + 1)) * x_bar_nk + (1 / (k + 1)) * x_k;
            t_bar_nk = (k / (k + 1)) * t_bar_nk + (1 / (k + 1)) * t_k;
            p_bar_nk = (k / (k + 1)) * p_bar_nk + (1 / (k + 1)) * p_k;
            y_bar_nk = (k / (k + 1)) * y_bar_nk + (1 / (k + 1)) * y_k;
            
            % -------------------- Calculate and Store Metrics --------------------
            % Calculate objective value using averaged p_bar_nk for consistency with dual methods
            p_over_v = p_bar_nk ./ v;
            p_over_v(p_over_v <= 0) = 1e-12; % Protect against non-positive values
            
            beta_i_inv_current = min(p_over_v, [], 2);
            beta_i_inv_current(beta_i_inv_current <= 0) = 1e-12; % Ensure positivity for log

            current_dual_obj = sum(p_bar_nk) - sum(B .* log(beta_i_inv_current));
            obj_values(iter_idx) = abs(current_dual_obj - fval_solver); % Gap to optimal dual value

            % Calculate distance to optimal prices using averaged p_bar_nk
            dis_p(iter_idx) = norm(p_bar_nk - p_opt_solver);

            % Check overall convergence (optional, usually for the outer loop)
            %%% Todo: we do not need to check the distance of dis_p
            if obj_values(iter_idx)  < epsilon % && dis_p(iter_idx) < epsilon
                break; % Early exit if converged
            end
        end % End inner loop

        % If inner loop broke early, total iterations will be less
        %%% Todo: we do not need to check the distance of dis_p
        if obj_values(iter_idx) < epsilon % && dis_p(iter_idx) < epsilon
            break;
        end

        % Restart the inner loop (Algorithm 1, Line 6)
        x_n0 = x_bar_nk;
        t_n0 = t_bar_nk;
        p_n0 = p_bar_nk;
        y_n0 = y_bar_nk;
    end % End outer loop

    time = toc; % End overall timer

    % Final solutions
    solution_x = x_bar_nk; % Return the averaged x
    solution_t = t_bar_nk; % Return the averaged t
    solution_p = p_bar_nk; % Return the averaged p
    solution_y = y_bar_nk; % Return the averaged y

    % Trim results arrays
    obj_values = obj_values(1:iter_idx);
    dis_p = dis_p(1:iter_idx);
    iter_total = iter_idx;

    % Plotting if plot_flag is true
    if plot_flag
        figure;
        iterations_vec = 0:(iter_total - 1); % Now includes 0 for initial point
        
        subplot(2, 1, 1);
        semilogy(iterations_vec, obj_values, 'LineWidth', 2);
        xlabel('Total Inner Iteration');
        ylabel('Dual Objective Value Gap');
        title('PDHG (Primal-Dual) - Objective Gap Convergence');
        grid on;

        subplot(2, 1, 2);
        semilogy(iterations_vec, dis_p, 'LineWidth', 2);
        xlabel('Total Inner Iteration');
        ylabel('Distance to P*');
        title('PDHG (Primal-Dual) - Price Distance Convergence');
        grid on;
    end
end