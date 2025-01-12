function [gap, gap_array, matrix_backup]= linear_compute_gap(v, B, mu)
    % Input:
    % v - parameter matrix v \in R^{n*m}
    % B - vector B \in R^{1*m}
    % mu - vector mu \in R^{1*m}
    % Output:
    % gap - scalar gap = min_{i} gap_i

    % Get the dimensions of the matrix
    % [n, m] = size(v);

    % Compute the index set and beyond
    matrix = log(v) - mu;
    matrix_backup = matrix;
    top_value = max(matrix,[],2);
    top_index = abs(matrix - top_value)<1e-5;
    matrix(top_index) = -inf;
    second_top_value = max(matrix,[],2);
    gap_array = top_value - second_top_value;
    gap = min(gap_array);
    % % Compute the gap for each i
    % for i = 1:n
    %     % Compute the log(v_{ij}) - mu_{j} for each j
    %     log_v_minus_mu = log(v(i, :)) - mu;
        
    %     % Find the maximum value and all indices where the maximum value occurs
    %     max_val = max(log_v_minus_mu);
    %     max_indices = log_v_minus_mu == max_val;
    %     % Set the maximum value to -inf for all indices where the maximum value occurs
    %     log_v_minus_mu(max_indices) = -inf;

        
    %     % Find the second maximum value
    %     second_max_val = max(log_v_minus_mu);
    %     % Compute the gap for this i
    %     gap_array(i) = max_val - second_max_val;
    % end
    % % Compute the minimum gap
    % gap = min(gap_array);
end