function solution = linear_class_generation(A, v, B)
    %%% Todo: consider the problem brought by numerical precision and log(v) 
    %%% Todo: consider change back to binary for classification
    %%% * Given binary matrix - output the solution
    % A: Input matrix with n rows and m columns, entries are either 0 or 1
    % v: Constant matrix with the same size as A
    % B: n by 1 vector
    % solution: 1*m vector

    [n, m] = size(A);
    
    % Step 1: Identify row classes and corresponding column classes
    visited_rows = false(1, n);
    solution = zeros(1, m);
    
    for i = 1:n
        if ~visited_rows(i)
            % Find all rows that share at least one 1 with row i
            same_rows = find(any(A & A(i, :), 2));
            visited_rows(same_rows) = true;
            
            % Initialize the column class for this row class
            col_class = zeros(length(same_rows), m);
            
            % Fill the column class matrix using logical indexing
            for j = 1:length(same_rows)
                row_idx = same_rows(j);
                col_class(j, A(row_idx, :) == 1) = log(v(row_idx, A(row_idx, :) == 1));
            end
            
            % Step 2: Calculate gap vector for this class
            [gap_vector, mu_basis] = calculate_gap_vector(col_class, B(same_rows));
            % Step 3: Calculate solution vector
            active_cols = find(gap_vector ~= -inf);
            solution(active_cols) = mu_basis + gap_vector(active_cols);
        end
    end
end

function [gap_vector, mu_basis] = calculate_gap_vector(matrix, B_subset)
    % matrix: Input matrix (d by m), entries are either 0 or nonzero
    % B_subset: Subset of B corresponding to the rows of the matrix
    % gap_vector: 1*m gap vector
    % mu_basis: Basis column value for the solution vector

    [d, m] = size(matrix);
    
    % Step 1: Choose basis and construct checked set
    active_cols = find(any(matrix, 1));
    basis_col = min(active_cols);
    checked_set = basis_col;
    gap_vector = -inf(1, m);
    gap_vector(basis_col) = 0;
    
    % Step 2: Start from the row that basis column is active in
    active_rows = find(matrix(:, basis_col) ~= 0);
    for i = 1:length(active_rows)
        row_idx = active_rows(i);
        other_active_cols = find(matrix(row_idx, :) ~= 0);
        for j = 1:length(other_active_cols)
            col_idx = other_active_cols(j);
            if ~ismember(col_idx, checked_set)
                gap_vector(col_idx) = matrix(row_idx, col_idx) - matrix(row_idx, basis_col);
                checked_set = [checked_set, col_idx];
            end
        end
    end
    
    % Step 3: Outer iteration loop
    while length(checked_set) < length(active_cols)
        % Find an activation index not in the checked set with gap not zero and not negative infinity
        candidate_cols = setdiff(active_cols, checked_set);
        candidate_cols = candidate_cols(gap_vector(candidate_cols) ~= 0 & gap_vector(candidate_cols) ~= -inf);
        if isempty(candidate_cols)
            break;
        end
        next_col = candidate_cols(1);
        
        % Inner iteration loop
        active_rows = find(matrix(:, next_col) ~= 0);
        for i = 1:length(active_rows)
            row_idx = active_rows(i);
            other_active_cols = find(matrix(row_idx, :) ~= 0);
            for j = 1:length(other_active_cols)
                col_idx = other_active_cols(j);
                if ~ismember(col_idx, checked_set)
                    gap_vector(col_idx) = gap_vector(next_col) + matrix(row_idx, col_idx) - matrix(row_idx, next_col);
                    checked_set = [checked_set, col_idx];
                end
            end
        end
    end
    
    % Calculate mu_basis
    mu_basis = log(sum(B_subset)) - log(sum(exp(gap_vector(gap_vector ~= -inf))));
end