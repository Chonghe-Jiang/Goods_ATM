function [solution, time] = linear_primal_yalmip(v, B)
    % Input:
    % v - parameter matrix v \in R^{n*m}
    % B - parameter vector B \in R^{n*1}
    % Output:
    % solution - solution of the primal EG problem
    % time - time taken to solve the problem

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Start timing
    tic;

    % Define the decision variable using YALMIP
    x = sdpvar(n, m, 'full');

    % Define the objective function
    objective = sum(B .* log(sum(v .* x, 2)));

    % Define the constraints
    constraints = [sum(x, 1) <= 1, x >= 0];

    % Set the solver to Gurobi (or any other solver supported by YALMIP)
    options = sdpsettings('solver', 'mosek');

    % Solve the problem
    diagnostics = optimize(constraints, -objective, options);

    % End timing
    time = toc;

    % Check if the optimization was successful
    if diagnostics.problem == 0
        % Extract the solution
        solution = value(x);
    else
        error('Optimization failed: %s', diagnostics.info);
    end
end