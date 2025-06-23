function gap = linear_optimal_gap(v, mu)
    matrix = log(v) - mu;
    top_value = max(matrix, [], 2);
    %%% ! Important step: machine accuracy
    matrix_backup = abs(matrix - top_value) < 1e-8; %%% Todo: need to revise - for numerical accuracy
    matrix(matrix_backup) = -inf;
    second_top_value = max(matrix, [], 2);
    gap_array = top_value - second_top_value;
    gap = min(gap_array);
end