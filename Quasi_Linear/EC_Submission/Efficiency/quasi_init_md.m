function b = quasi_init_md(p, B)
    % Make sure that we satisfy the case: \sum p <= \sum B
    n = size(B, 1);
    m = size(p, 2);

    % Constraint 1: sum_i b_{i j} = p_j (equality constraint)
    Aeq = kron(eye(m), ones(1, n)); % sum_i b_{i j} = p_j
    beq = p';
    
    % Inequality constraint: sum_j b_{i j} <= B_i for all i
    A = kron(ones(1, m), eye(n)); % sum_j b_{i j} <= B_i
    b_ineq = B; % B_i

    % Objective function (can be set arbitrarily)
    f = zeros(n * m, 1);
    
    % Variable bounds
    lb = zeros(n * m, 1);
    ub = inf(n * m, 1);

    % Solve the linear programming problem
    b_flat = linprog(f, A, b_ineq, Aeq, beq, lb, ub);

    % Reshape the result back to an n x m matrix
    b = reshape(b_flat, n, m);
end
