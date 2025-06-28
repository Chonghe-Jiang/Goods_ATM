function [solution_x, solution_t, solution_p, solution_y, obj_values, dis_p, time, iter_total] = linear_primal_dual_pdhg_adaptive(v, B, x0, t0, p0, y0, max_outer_iter, K_inner_iter, eta_initial, omega_initial, theta, epsilon, plot_flag, p_opt_solver, fval_solver)
    % linear_primal_dual_pdhg_adaptive - Implements PDHG with adaptive step-sizes for eta and omega.
    % Based on the PDHG algorithm in "PDHCG: A Scalable First-Order Method..." (arXiv:2506.06258v1).
    % This version incorporates an adaptive step-size strategy for sigma and tau.

    % Input:
    % v              - Valuation matrix u_ij (n x m)
    % B              - Budget vector w_i (n x 1)
    % x0             - Initial allocation matrix x (n x m)
    % t0             - Initial auxiliary variable t (n x 1)
    % p0             - Initial price vector p (1 x m)
    % y0             - Initial dual variable y (n x 1)
    % max_outer_iter - Maximum number of outer loop iterations (restarts)
    % K_inner_iter   - Number of inner iterations per outer loop
    % eta_initial    - Initial overall step-size
    % omega_initial  - Initial primal weight (balances primal/dual steps)
    % theta          - Adaptation factor for omega and eta
    % epsilon        - Tolerance for convergence
    % plot_flag      - True to plot, false to not plot
    % p_opt_solver   - Optimal price solution from an external solver (for plotting)
    % fval_solver    - Optimal primal objective value from an external solver (for plotting)

    [n, m] = size(v);

    % --- Adaptive Step-Size Parameters ---
    eta = eta_initial;
    omega = omega_initial;
    % Define bounds for primal weight to prevent extreme values
    omega_min = 0.1;
    omega_max = 10;
    
    % Initialize variables for outer loop
    x_n0 = x0;
    t_n0 = t0;
    p_n0 = p0;
    y_n0 = y0;

    % Store results
    obj_values = zeros(max_outer_iter * K_inner_iter + 1, 1);
    dis_p = zeros(max_outer_iter * K_inner_iter + 1, 1);
    iter_idx = 1;

    tic; % Start overall timer

    % --- Record Initial Objective Function Value and Distance ---
    p_over_v_initial = p_n0 ./ v;
    p_over_v_initial(p_over_v_initial <= 0) = 1e-12;
    beta_i_inv_initial = min(p_over_v_initial, [], 2);
    beta_i_inv_initial(beta_i_inv_initial <= 0) = 1e-12;
    initial_dual_obj = sum(p_n0) - sum(B .* log(beta_i_inv_initial));
    obj_values(iter_idx) = abs(initial_dual_obj - fval_solver);
    dis_p(iter_idx) = norm(p_n0 - p_opt_solver);
    % -----------------------------------------------------------

    for outer_n = 0:max_outer_iter-1 % Outer loop: Restart mechanism
        
        % --- Calculate adaptive step-sizes for this outer loop ---
        tau_step = eta / omega;
        sigma_step = eta * omega;
        
        % Initialize inner loop state
        x_bar_nk = x_n0;
        t_bar_nk = t_n0;
        p_bar_nk = p_n0;
        y_bar_nk = y_n0;

        % Keep track of previous iterates for extrapolation
        x_prev = x_n0;
        t_prev = t_n0;
        p_prev = p_n0;
        y_prev = y_n0;

        % Set current iterates for the first inner loop update
        x_k = x_n0;
        t_k = t_n0;
        p_k = p_n0;
        y_k = y_n0;

        for k = 0:K_inner_iter-1 % Inner loop
            iter_idx = iter_idx + 1;

            % Extrapolated primal variables
            x_extrapolated = 2 * x_k - x_prev;
            t_extrapolated = 2 * t_k - t_prev;

            % --- Dual Variables Update (p, y) ---
            p_new = p_k + sigma_step * (sum(x_extrapolated, 1) - 1);
            y_new = y_k + sigma_step * (t_extrapolated - sum(v .* x_extrapolated, 2));

            % --- Primal Variables Update (x, t) ---
            x_new_unprojected = x_k - tau_step * (p_new - v .* y_new);
            x_new = max(x_new_unprojected, 0);

            sqrt_term_t = (tau_step * y_new - t_k).^2 + 4 * tau_step * B;
            t_new = (t_k - tau_step * y_new + sqrt(sqrt_term_t)) / 2;
            t_new = max(t_new, 1e-12);

            % Store current and previous iterates for next iteration
            x_prev = x_k; t_prev = t_k; p_prev = p_k; y_prev = y_k;
            x_k = x_new;  t_k = t_new;  p_k = p_new;  y_k = y_new;

            % --- Update Averaged Solutions ---
            x_bar_nk = (k / (k + 1)) * x_bar_nk + (1 / (k + 1)) * x_k;
            t_bar_nk = (k / (k + 1)) * t_bar_nk + (1 / (k + 1)) * t_k;
            p_bar_nk = (k / (k + 1)) * p_bar_nk + (1 / (k + 1)) * p_k;
            y_bar_nk = (k / (k + 1)) * y_bar_nk + (1 / (k + 1)) * y_k;
            
            % --- Calculate and Store Metrics ---
            p_over_v = p_bar_nk ./ v;
            p_over_v(p_over_v <= 0) = 1e-12;
            beta_i_inv_current = min(p_over_v, [], 2);
            beta_i_inv_current(beta_i_inv_current <= 0) = 1e-12;
            current_dual_obj = sum(p_bar_nk) - sum(B .* log(beta_i_inv_current));
            obj_values(iter_idx) = abs(current_dual_obj - fval_solver);
            dis_p(iter_idx) = norm(p_bar_nk - p_opt_solver);

            if obj_values(iter_idx) < epsilon
                break;
            end
        end % End inner loop

        % --- Adaptive Step-Size and Primal Weight Update ---
        % This check is performed after each restart cycle.
        if mod(outer_n, 1) == 0 && outer_n > 0
            % Calculate primal residual (constraint violation)
            primal_res1 = norm(sum(x_bar_nk, 1) - 1, inf);
            primal_res2 = norm(t_bar_nk - sum(v .* x_bar_nk, 2), inf);
            r_primal = max(primal_res1, primal_res2);

            % Calculate dual residual
            grad_L_x = p_bar_nk - v .* y_bar_nk;
            r_dual = norm(x_bar_nk .* grad_L_x, 1) / (1 + norm(x_bar_nk,1));
            
            % Adjust omega to balance primal and dual progress
            if r_primal > (1 + theta) * r_dual
                omega = omega / (1 - theta); % Primal residual is larger, increase primal step (decrease omega)
            elseif r_dual > (1 + theta) * r_primal
                omega = omega * (1 - theta); % Dual residual is larger, increase dual step (increase omega)
            end
            % Enforce bounds on omega
            % omega = max(min(omega, omega_max), omega_min);

            % --- NEW: Adjust and Bound ETA ---
            % ! Check the operation here, whether it is necessary to adjust eta based on residuals.
            % This logic was added based on the paper's description to control the overall step size.
            % if r_primal > (1 + theta) * r_dual || r_dual > (1 + theta) * r_primal
            %     % If residuals are imbalanced, be more conservative with the overall step size.
            %     eta = eta * (1 - theta / 2);
            % else
            %     % If residuals are balanced, we can be more aggressive.
            %     eta = eta / (1 - theta / 2);
            % end
            % Apply the bounds as specified in the paper
            eta = max(0.01 * eta_initial, min(eta, 3 * eta_initial));
            % --- END NEW LOGIC ---
        end


        if obj_values(iter_idx) < epsilon
            break;
        end

        % Restart the inner loop
        x_n0 = x_bar_nk;
        t_n0 = t_bar_nk;
        p_n0 = p_bar_nk;
        y_n0 = y_bar_nk;
    end % End outer loop

    time = toc;

    % Final solutions
    solution_x = x_bar_nk;
    solution_t = t_bar_nk;
    solution_p = p_bar_nk;
    solution_y = y_bar_nk;

    % Trim results arrays
    obj_values = obj_values(1:iter_idx);
    dis_p = dis_p(1:iter_idx);
    iter_total = iter_idx;

    % Plotting if flag is true
    if plot_flag
        figure;
        iterations_vec = 0:(iter_total - 1);
        
        subplot(2, 1, 1);
        semilogy(iterations_vec, obj_values, 'LineWidth', 2);
        xlabel('Total Inner Iteration');
        ylabel('Dual Objective Value Gap');
        title('PDHG (Adaptive) - Objective Gap Convergence');
        grid on;

        subplot(2, 1, 2);
        semilogy(iterations_vec, dis_p, 'LineWidth', 2);
        xlabel('Total Inner Iteration');
        ylabel('Distance to P*');
        title('PDHG (Adaptive) - Price Distance Convergence');
        grid on;
    end
end
