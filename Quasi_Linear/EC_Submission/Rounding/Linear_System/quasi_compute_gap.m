function [gap, gap_array, activation_matrix] = quasi_compute_gap(v, B, mu, radius_current)
    [n, m] = size(v);
    matrix = log(v) - mu;
    top_value = max(matrix, [], 2);
    %%% ! Important step: machine accuracy
    matrix_backup = abs(matrix - top_value) < 1e-4; %%% Todo: need to revise - for numerical accuracy
    matrix(matrix_backup) = -inf;
    second_top_value = max(matrix, [], 2);
    gap_array = top_value - second_top_value;
    gap = min(gap_array);
    delta_threshold = radius_current*2;
    matrix = log(v) - mu;
    activation_matrix = (matrix - top_value) >= -delta_threshold;
end