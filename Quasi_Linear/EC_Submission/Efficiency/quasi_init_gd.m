function p = quasi_init_gd(l, u, c)
    %%% ! High level -> determine a value and then make md to satisfy this
    % Input:
    % l: lower bound vector (1 x m)
    % u: upper bound vector (1 x m)
    % c: constant sum of vector elements
    % Output:
    % x: generated vector (1 x m)

    m = length(l);
    
    % Construct the linear programming constraint matrix and vector
    A = ones(1, m); % Corresponds to sum_i x_i <= c
    b = c;
    
    % Linear programming objective function (can be set arbitrarily)
    f = zeros(m, 1);
    
    % Linear programming variable lower and upper bounds
    lb = l';
    ub = u';
    %%% ! Change it to quasi linear case with inequality
    % Use linprog to solve the linear programming problem
    p = linprog(f, A, b, [], [], lb, ub);
    
    % Reshape the result back to a 1 x m vector
    p = p';
end