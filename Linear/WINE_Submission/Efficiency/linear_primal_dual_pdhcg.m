function [solution_x, solution_p, obj_values, dis_p, time, iter_total] = linear_primal_dual_pdhcg(v, B, x0, p0, max_outer_iter, K_inner_iter, sigma_step, tau_step, epsilon, plot_flag, p_opt_solver, fval_solver)
    % linear_primal_dual_pdhcg - Implements Primal-Dual Hybrid Conjugate Gradient (PDHCG) for Fisher Market Equilibrium.
    % Based on Algorithm 2 and equations (5), (10), (11) from "PDHCG: A Scalable First-Order Method..." (arXiv:2506.06258v1).
    % This version uses bisection search for the inner x-subproblem.

    % Input:
    % v              - Valuation matrix u_ij (n x m)
    % B              - Budget vector w_i (n x 1)
    % x0             - Initial allocation matrix x (n x m)
    % p0             - Initial price vector p (1 x m)
    % max_outer_iter - Maximum number of outer loop iterations (restarts)
    % K_inner_iter   - Number of inner iterations per outer loop
    % sigma_step     - Dual step size (sigma)
    % tau_step       - Primal step size (tau)
    % epsilon        - Tolerance for convergence (based on primal objective gap)
    % plot_flag      - True to plot, false to not plot
    % p_opt_solver   - Optimal price solution from external solver
    % fval_solver    - Optimal primal objective value from external solver

    % Output:
    % solution_x     - Final allocation matrix x
    % solution_p     - Final price vector p
    % obj_values     - Array of primal objective function values (gap to optimal)
    % dis_p          - Array of distances between current p and p_opt_solver
    % time           - Time taken to solve the problem
    % iter_total     - Total number of inner iterations performed

    [n, m] = size(v);

    % Initialize variables for outer loop (Algorithm 2)
    x_n0 = x0; % x^(n,0)
    p_n0 = p0; % p^(n,0)
    
    % Ensure initial x is non-negative and sums to 1 per column.
    % This part is commented out because it's assumed x0 comes from linear_init_md
    % which handles initial feasibility. However, robust code might include it.
    % for j = 1:m
    %     if sum(x_n0(:,j)) > 0
    %         x_n0(:,j) = x_n0(:,j) / sum(x_n0(:,j));
    %     else
    %         x_n0(:,j) = ones(n, 1) / n; % Default to equal allocation if column sum is zero
    %     end
    %     x_n0(:,j) = max(x_n0(:,j), 1e-12); % Ensure strict positivity for log
    % end
    
    % Ensure p0 is a row vector and strictly positive.
    % This part is commented out because p0 from linear_init_gd is assumed feasible.
    % if iscolumn(p_n0)
    %     p_n0 = p_n0';
    % end
    % p_n0 = max(p_n0, 1e-12); % Ensure positivity

    % Store results (objective gap and price distance)
    obj_values = zeros(max_outer_iter * K_inner_iter, 1);
    dis_p = zeros(max_outer_iter * K_inner_iter, 1);
    iter_idx = 0;

    tic; % Start overall timer

    for outer_n = 0:max_outer_iter-1 % Outer loop: Restart mechanism
        % Initialize inner loop (Algorithm 2, Line 2)
        x_bar_nk = x_n0; % x_bar^(n,0)
        p_bar_nk = p_n0; % p_bar^(n,0)

        % Keep track of previous iterates for (2x^k - x^{k-1}) term
        % For k=0, use x_n0, p_n0 as initial values for extrapolation.
        % This is a common practice when x^{k-1} is not explicitly available.
        x_prev = x_n0;
        p_prev = p_n0;

        for k = 0:K_inner_iter-1 % Inner loop
            iter_idx = iter_idx + 1; % Total inner iteration counter

            % -------------------- Dual Variable Update (p) --------------------
            % (Algorithm 2, Line 3) and Eq. (5) for p.
            % Let's denote current iterate as x_k, p_k. Previous as x_prev, p_prev.
            % For the first inner iteration (k=0), x_k and p_k are x_n0 and p_n0.
            % For subsequent inner iterations (k>0), x_k and p_k are the results from the previous step.
            if k == 0
                x_k = x_n0;
                p_k = p_n0;
            end
            
            x_extrapolated = 2 * x_k - x_prev; % (2x^k - x^{k-1}) term
            
            % Update p (Equation 5)
            % sum_{i \in N} x_{ij} is sum(x_extrapolated) over i for each j.
            p_new = p_k + sigma_step * (sum(x_extrapolated, 1) - 1); % sum(x,1) sums columns
            p_new = max(p_new, 1e-12); % Ensure non-negativity and avoid log issues

            % -------------------- Primal Variable Update (x) --------------------
            % (Algorithm 2, Line 4) and Algorithm 3 (Bisection Search) for each row of x.
            x_new = zeros(n, m);
            for i = 1:n % Parallel update each row of x
                % Solve the subproblem for x_i using bisection search
                x_new(i, :) = solve_x_subproblem_bisection(v(i, :), B(i), p_new, x_k(i, :), tau_step);
            end
            x_new = max(x_new, 1e-12); % Ensure strict positivity for log in next iter

            % Store current and previous iterates for next iteration's extrapolation
            x_prev = x_k;
            p_prev = p_k;
            x_k = x_new;
            p_k = p_new;

            % -------------------- Update Averaged Solutions --------------------
            % (Algorithm 2, Line 5)
            x_bar_nk = (k / (k + 1)) * x_bar_nk + (1 / (k + 1)) * x_k;
            p_bar_nk = (k / (k + 1)) * p_bar_nk + (1 / (k + 1)) * p_k;
            
            % -------------------- Calculate and Store Metrics --------------------
            % Calculate objective value using averaged p_bar_nk for consistency with dual methods
            % Dual objective: sum(p) - sum(B .* log(min(p ./ v, [], 2)))
            p_over_v = p_bar_nk ./ v;
            p_over_v(p_over_v <= 0) = 1e-12; % Protect against non-positive values
            
            beta_i_inv_current = min(p_over_v, [], 2);
            % ! Check points
            % beta_i_inv_current(beta_i_inv_current <= 0) = 1e-12; % Ensure positivity for log

            current_dual_obj = sum(p_bar_nk) - sum(B .* log(beta_i_inv_current));
            obj_values(iter_idx) = abs(current_dual_obj - fval_solver); % Gap to optimal dual value (which equals primal optimal due to strong duality)

            % Calculate distance to optimal prices using averaged p_bar_nk
            dis_p(iter_idx) = norm(p_bar_nk - p_opt_solver);

            % Check overall convergence (optional, usually for the outer loop)
            if obj_values(iter_idx) < epsilon && dis_p(iter_idx) < epsilon
                break; % Early exit if converged
            end
        end % End inner loop

        % If inner loop broke early, total iterations will be less
        if obj_values(iter_idx) < epsilon && dis_p(iter_idx) < epsilon
            break;
        end

        % Restart the inner loop (Algorithm 2, Line 6)
        x_n0 = x_bar_nk;
        p_n0 = p_bar_nk;
    end % End outer loop

    time = toc; % End overall timer

    % Final solutions
    solution_x = x_bar_nk; % Return the averaged x
    solution_p = p_bar_nk; % Return the averaged p

    % Trim results arrays
    obj_values = obj_values(1:iter_idx);
    dis_p = dis_p(1:iter_idx);
    iter_total = iter_idx;

    % Plotting if plot_flag is true
    if plot_flag
        figure;
        iterations_vec = 1:iter_total;
        
        subplot(2, 1, 1);
        semilogy(iterations_vec, obj_values, 'LineWidth', 2);
        xlabel('Total Inner Iteration');
        ylabel('Dual Objective Value Gap'); % Changed Y-label
        title('PDHCG (Primal-Dual) - Objective Gap Convergence');
        grid on;

        subplot(2, 1, 2);
        semilogy(iterations_vec, dis_p, 'LineWidth', 2);
        xlabel('Total Inner Iteration');
        ylabel('Distance to P*');
        title('PDHCG (Primal-Dual) - Price Distance Convergence');
        grid on;
    end
end

% Helper function to solve x_i subproblem using bisection search (Algorithm 3)
% Solves for x_i (1 x m row vector)
function x_i_new = solve_x_subproblem_bisection(u_i, w_i, p_k_plus_1, x_i_k, tau_k)
    % u_i (1 x m) - i-th row of valuation matrix
    % w_i - i-th buyer's budget
    % p_k_plus_1 (1 x m) - current prices p^(k+1)
    % x_i_k (1 x m) - previous allocation x_i^k
    % tau_k - primal step size

    s_tol = 1e-8; % Tolerance for s (utility) in bisection

    % Algorithm 3, Line 1-3
    s_trial_init = sum(u_i .* x_i_k); % Initial guess for s (utility)

    % Compute x_tilde from s_trial_init using Eq. (11)
    % x_ij = proj_R+(x_ij^k - tau_k * p_j^{k+1} + tau_k * w_i / s)
    s_trial_init = max(s_trial_init, 1e-12); % Avoid division by zero or log(0)
    term_s = tau_k * w_i / s_trial_init;
    x_tilde_init_unproj = x_i_k - tau_k * p_k_plus_1 + term_s;
    x_tilde_init = max(x_tilde_init_unproj, 0); % Projection onto R+

    s_tilde_init = sum(u_i .* x_tilde_init); % s_tilde = u_i^T * x_tilde
    s_tilde_init = max(s_tilde_init, 1e-12);

    L = min(s_trial_init, s_tilde_init); % Lower bound for s*
    U = max(s_trial_init, s_tilde_init); % Upper bound for s*

    % If L or U are zero or very small, adjust them
    if U <= 1e-12
        U = 1; % Default to a reasonable upper bound for utility
        L = 1e-12; % Keep L positive
    end
    if L > U % Should not happen with min/max, but safety check
        temp = L; L = U; U = temp;
    end

    % Algorithm 3, Line 4-7: Bisection search for s*
    max_bisection_iter = 100; % Max iterations for inner bisection loop
    bisection_iter = 0;
    
    while (U - L > s_tol) && (bisection_iter < max_bisection_iter)
        bisection_iter = bisection_iter + 1;
        
        s_mid = (L + U) / 2; % Algorithm 3, Line 5
        s_mid = max(s_mid, 1e-12); % Ensure positivity

        % Compute x_tilde from s_mid using Eq. (11)
        term_s_mid = tau_k * w_i / s_mid;
        x_tilde_mid_unproj = x_i_k - tau_k * p_k_plus_1 + term_s_mid;
        x_tilde_mid = max(x_tilde_mid_unproj, 0); % Projection onto R+

        s_tilde_mid = sum(u_i .* x_tilde_mid); % s_tilde = u_i^T * x_tilde_mid
        s_tilde_mid = max(s_tilde_mid, 1e-12);

        % Update L and U based on (s_mid, s_tilde_mid) pair (Algorithm 3, Line 6)
        % This is the "enhanced bisection" from Remark 1
        U = min(U, max(s_mid, s_tilde_mid));
        L = max(L, min(s_mid, s_tilde_mid));
    end

    % After bisection, U and L are close to s*. We can use (L+U)/2 as s*.
    current_s_star = (L + U) / 2;
    current_s_star = max(current_s_star, 1e-12); % Final check for positivity

    % Compute the final x_i_new using the found s* (Algorithm 3, Output)
    term_s_star = tau_k * w_i / current_s_star;
    x_i_new_unproj = x_i_k - tau_k * p_k_plus_1 + term_s_star;
    x_i_new = max(x_i_new_unproj, 0); % Final projection onto R+
    
    % Final check: Although not explicitly part of Algo 3's output, it's good practice
    % to ensure the sum(u_i .* x_i_new) is close to current_s_star.
    % The primal subproblem aims to find x_i_new such that u_i' * x_i_new = s*.
    % The column sum constraint (sum_i x_ij = 1) is implicitly handled by the dual variable p_j's update.
end