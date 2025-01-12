function b = linear_init_md(p, B)
    % Make sure that we satisfy the case: \sum p <= \sum B
    n = size(B, 1);
    m = size(p, 2);

    % Constraint 1: sum_i b_{i j} = p_j (equality constraint)
    Aeq1 = kron(eye(m), ones(1, n)); % sum_i b_{i j} = p_j
    beq1 = p';
    
    % Constraint 2: sum_j b_{i j} = B_i (equality constraint)
    Aeq2 = repmat(eye(n), 1, m); % sum_j b_{i j} = B_i
    beq2 = B;
    
    % Combine equality constraints
    Aeq = [Aeq1; Aeq2];
    beq = [beq1; beq2];

    % Inequality constraint: sum p <= sum B
    A = ones(1, n * m); % sum_{i, j} b_{i j} <= sum B
    b = sum(B); % sum B

    % Objective function (can be set arbitrarily)
    f = zeros(n * m, 1);
    
    % Variable bounds
    lb = zeros(n * m, 1);
    ub = inf(n * m, 1);

    % Solve the linear programming problem
    b_flat = linprog(f, A, b, Aeq, beq, lb, ub);

    % Reshape the result back to an n x m matrix
    b = reshape(b_flat, n, m);
end