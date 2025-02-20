function [oracle_result, optimal_ce, optimal_value, sum_exp_mu] = linear_exact_oracle(v, B, mu, radius_current)
    % linear_exact_oracle: Computes the exact oracle result and optimal certificate
    % Inputs:
    %   v: n x m matrix
    %   B: n x 1 vector
    %   mu: 1 x m vector
    % Outputs:
    %   oracle_result: True if the solution is feasible, False otherwise
    %   optimal_ce: Optimal certificate (solution) if feasible, NaN otherwise

    %%% ! Step 1: Compute the gap and activation matrix
    [~ ,~ , activation_matrix] = linear_compute_gap(v, B, mu, radius_current);

    %%% ! Step 2: Generate the solution using linear_class_generation: two stages inside
    solution = linear_class_generation(activation_matrix, v, B);
    solution
    % Step 3: Check feasibility using linear_max_flow
    [is_feasible, optimal_value, sum_exp_mu] = linear_max_flow(solution, B, v, activation_matrix);
    % Step 4: Return the appropriate results based on feasibility
    if is_feasible
        oracle_result = true;
        optimal_ce = solution;
    else
        oracle_result = false;
        optimal_ce = NaN;
    end
end