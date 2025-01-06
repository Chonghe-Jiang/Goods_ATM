function gap= linear_compute_gap_cheating(v, B, mu)
    %%% ! This is to compute the \Delta in computation, also output a binary matrix matrix_backup
    %%% ! In addition, you also need to output 
    % Input:
    % v - parameter matrix v \in R^{n*m}
    % B - vector B \in R^{1*m}
    % mu - vector mu \in R^{1*m}
    % Output:
    % gap - scalar gap = min_{i} gap_i

    % Get the dimensions of the matrix
    [n, m] = size(v);

    % Compute the index set and beyond
    matrix = log(v) - mu;
    
    % Find the largest element in each row
    top_value = max(matrix, [], 2);
    
    % Create a binary matrix where the largest element in each row is 1
    % Use a numerical precision threshold of 1e-8
    matrix_backup = abs(matrix - top_value) < 1e-4; %%% Todo: need to revise
    
    % Set the largest elements to -inf to find the second largest
    matrix(matrix_backup) = -inf;
    
    % Find the second largest element in each row
    second_top_value = max(matrix, [], 2);
    
    % Compute the gap array
    gap_array = top_value - second_top_value;
    
    % Compute the minimum gap
    gap = min(gap_array);
end