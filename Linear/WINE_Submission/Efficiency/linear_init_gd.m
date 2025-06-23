function p = linear_init_gd(l, u, c)
    % Input:
    % l: lower bound vector (1 x m)
    % u: upper bound vector (1 x m)
    % c: constant sum of vector elements
    % Output:
    % x: generated vector (1 x m)

    m = length(l);
    
    % Construct the linear programming constraint matrix and vector
    Aeq = ones(1, m); % Corresponds to sum_i x_i = c
    beq = c;
    
    % Linear programming objective function (can be set arbitrarily)
    f = zeros(m, 1);
    
    % Linear programming variable lower and upper bounds
    lb = l';
    ub = u';
    
    % Use linprog to solve the linear programming problem
    p = linprog(f, [], [], Aeq, beq, lb, ub);
    
    % Reshape the result back to a 1 x n vector
    p = p';
end
