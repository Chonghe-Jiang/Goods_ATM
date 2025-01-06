function [solution, time] =linear_primal_cvx(v, B)
    % Input:
    % v - parameter matrix v \in R^{n*m}
    % Output:
    % solution - solution of the primal EG problem
    % time - time taken to solve the problem

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Start timing
    tic;

    % Define the EG problem using CVX
    cvx_begin
        variable x(n, m)
        maximize(sum(B .* log(sum(v .* x, 2))))
        subject to
            sum(x, 1) <= 1;
            x >= 0;
    cvx_end

    % End timing
    time = toc;

    % Extract the solution
    solution = x;
end