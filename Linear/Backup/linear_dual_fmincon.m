function [p_opt, fval] = linear_dual_fmincon(v, B, p_0)
    % Inputs:
    % v - n x m matrix
    % B - n x 1 vector
    % p_0 - 1 x m vector (initial guess for p)
    
    % Outputs:
    % p_opt - 1 x m vector (optimal p)
    % fval - scalar (optimal objective function value)
    p_0 = p_0';
    [n, m] = size(v);
    
    % Find an initial beta_0 within the constraints
    beta_0 = find_initial_beta(v, p_0);
    
    % Combine p_0 and beta_0 into a single vector
    x0 = [p_0; beta_0];
    
    % Define the objective function
    fun = @(x) sum(x(1:m)) - sum(B .* log(x(m+1:end)));
    
    % Define the linear inequality constraint matrix A
    A = zeros(n*m, m+n);
    for i = 1:n
        for j = 1:m
            row = (i-1)*m + j;
            A(row, j) = -1;
            A(row, m+i) = v(i, j);
        end
    end
    
    % Define the nonlinear constraint function
    nonlcon = @(x) deal([], A*x);
    
    % Define the lower bounds for the variables
    lb = [-inf(m, 1); zeros(n, 1)];
    
    % Use fmincon to solve the optimization problem
    % options = optimoptions('fmincon', 'Display', 'iter');
    [x, fval] = fmincon(fun, x0, [], [], [], [], lb, [], nonlcon, options);
    
    % Extract the optimal p
    p_opt = x(1:m)';
end

function beta_0 = find_initial_beta(v, p_0)
    % This function finds an initial beta_0 within the constraints
    % p_j >= v_{ij} * beta_i
    
    [n, m] = size(v);
    beta_0 = zeros(n, 1);
    
    % For each i, find the minimum value of p_j / v_{ij} and set beta_0(i)
    % to a value that satisfies the constraints
    for i = 1:n
        valid_ratios = p_0 ./ v(i, :)';
        valid_ratios(v(i, :) == 0) = inf; % Avoid division by zero
        beta_0(i) = min(valid_ratios);
        if beta_0(i) <= 0
            beta_0(i) = 1; % Arbitrary positive value
        end
    end
end