function [oracle_result, optimal_ce, optimal_value, sum_exp_mu] = quasi_exact_oracle(v, B, mu, radius_current)
    % linear_exact_oracle: Computes the exact oracle result and optimal certificate
    % Inputs:
    %   v: n x m matrix
    %   B: n x 1 vector
    %   mu: 1 x m vector
    % Outputs:
    %   oracle_result: True if the solution is feasible, False otherwise
    %   optimal_ce: Optimal certificate (solution) if feasible, NaN otherwise

    %%% ! Step 1: Compute the gap and activation matrix
    %%% ! Step 2: Quasi-Linear Utility modification - we should input new v and new mu
    v = [ones(length(B),1),v]; % add a column
    mu = [0,mu]; % add an element
    [~ ,~ , activation_matrix] = quasi_compute_gap(v, B, mu, radius_current);

    %%% ! Step 2: Generate the solution using linear_class_generation: two stages inside
    solution = quasi_class_generation(activation_matrix, v, B);
    % Step 3: Check feasibility using linear_max_flow
    solution = [0,solution]; % Todo: check whether reasonable or not to do this for solution
    [is_feasible, optimal_value, sum_exp_mu] = quasi_max_flow(solution, B, v, activation_matrix);
    % Step 4: Return the appropriate results based on feasibility
    if is_feasible
        oracle_result = true;
        solution = solution(2:end); % remove the first element
        optimal_ce = solution;
    else
        oracle_result = false;
        optimal_ce = NaN;
    end
end