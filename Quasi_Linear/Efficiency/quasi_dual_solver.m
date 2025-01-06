function [p_opt, beta_opt, objective_value, solve_time] = quasi_dual_solver(n, m, B, v)
    %%% ! we change the problem formulation, including a new inequality constrain

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
        constraints = [constraints, beta(i) <= 1]; % Add beta_i <= 1 constraint
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