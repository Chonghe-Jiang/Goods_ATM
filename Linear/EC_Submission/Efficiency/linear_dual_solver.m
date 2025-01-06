function [p_opt, beta_opt, objective_value, solve_time] = linear_dual_solver(n, m, B, v)
    %LINEAR_DUAL_SOLVER Solves the optimization problem with given parameters
    %   [p_opt, beta_opt, objective_value, solve_time] = linear_dual_solver(n, m, B, v)
    %   Inputs:
    %       n - Range of i
    %       m - Range of j
    %       B - Values for B_i
    %       v - Values for v_{ij}
    %   Outputs:
    %       p_opt - Optimal values for p
    %       beta_opt - Optimal values for beta
    %       objective_value - Optimal objective function value
    %       solve_time - Time taken to solve the optimization problem

    % Define variables
    p = sdpvar(m, 1); % p_j variables
    beta = sdpvar(n, 1); % beta_i variables

    % Define objective function
    objective = sum(p) - sum(B .* log(beta));

    % Define constraints
    constraints = [];
    for i = 1:n
        for j = 1:m
            constraints = [constraints, p(j) >= v(i, j) * beta(i)];
        end
    end

    % Set optimization problem
    options = sdpsettings('solver', 'mosek', 'verbose', 0);

    % Record start time
    start_time = tic;

    % Solve optimization problem
    sol = optimize(constraints, objective, options);

    % Record end time
    end_time = toc(start_time);

    % Calculate solve time
    solve_time = end_time;

    % Check solution status
    if sol.problem == 0
        % Extract solutions
        p_opt = value(p);
        beta_opt = value(beta);
        objective_value = value(objective);
    else
        error('Optimization failed with error message: %s', sol.info);
    end
end