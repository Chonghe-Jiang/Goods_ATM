function b = quasi_init_md(p, B)
    % Make sure that we satisfy the case: \sum p = \sum B
    n = size(B, 1);
    m = size(p, 2);

    % Constraints for sum_i b_{i j} = p_j for j = 1:m
    Aeq1 = kron(eye(m), ones(1, n)); % sum_i b_{i j} = p_j
    beq1 = p';
    
    % Constraints for sum_j b_{i j} + δ_i = B_i
    Aeq2 = [repmat(eye(n), 1, m), eye(n)]; % sum_j b_{i j} + δ_i = B_i
    beq2 = B;
    
    Aeq = [Aeq1, zeros(size(Aeq1, 1), n); Aeq2];
    beq = [beq1; beq2];

    f = zeros(n*(m+1), 1);
    
    lb = zeros(n*(m+1), 1);
    ub = inf(n*(m+1), 1);

    b_flat = linprog(f, [], [], Aeq, beq, lb, ub);

    b = reshape(b_flat, n, m+1);
end
