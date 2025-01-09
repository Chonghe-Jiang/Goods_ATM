function solution = linear_class_generation(A, v, B)
    %%% * Given binary matrix - output the solution
    % A: Input matrix with n rows and m columns, entries are either 0 or 1
    % v: Constant matrix with the same size as A
    % B: n by 1 vector
    % solution: 1*m vector
    %%% ! Step 1: Do classification of different classes
    %%% Todo: the input A -> output classes; we use brand new function type to represent it 
    [row_classes, column_classes, col_class_matrices] = linear_return_activation(A, v);
    %%% ! Step 2: Do computation for every class - use the above information
    m = size(A, 2); % Number of columns in A
    solution = zeros(1, m); % Initialize solution vector

    for i = 1:length(row_classes)
        % Extract same_rows and col_class for the current class
        same_rows = row_classes{i}; % Rows in the current class
        col_class = col_class_matrices{i}; % Derived matrix for the current class
        B_class = B(same_rows, :);
        [gap_vector, mu_basis] = calculate_gap_vector(col_class, B_class);
        active_cols = find(gap_vector ~= -inf);
        solution(active_cols) = mu_basis + gap_vector(active_cols);
    end
    
    %%% Old Approach
    % [n, m] = size(A);
    % visited_rows = false(1, n);
    % solution = zeros(1, m);
    
    % for i = 1:n
    %     if ~visited_rows(i)
    %         % Find all rows that share at least one 1 with row i
    %         %%%  Operating to find rows
    %         same_rows = find(any(A & A(i, :), 2));
    %         visited_rows(same_rows) = true;
            
    %         %%%  Operating to find columns
    %         % Initialize the column class for this row class
    %         col_class = zeros(length(same_rows), m);
            
    %         % Fill the column class matrix using logical indexing
    %         for j = 1:length(same_rows)
    %             row_idx = same_rows(j);
    %             col_class(j, A(row_idx, :) == 1) = log(v(row_idx, A(row_idx, :) == 1));
    %         end
            
    %         % Step 2: Calculate gap vector for this class
    %         [gap_vector, mu_basis] = calculate_gap_vector(col_class, B(same_rows));
    %         % Step 3: Calculate solution vector
    %         active_cols = find(gap_vector ~= -inf);
    %         solution(active_cols) = mu_basis + gap_vector(active_cols);
    %     end
    % end
end
%%% ! This is a new approach for the gap calculation
function [gap_vector, mu_basis] = calculate_gap_vector(matrix, B_subset)
    % matrix: Input matrix (d by m), entries are either 0 or nonzero
    % B_subset: Subset of B corresponding to the rows of the matrix
    % gap_vector: 1*m gap vector
    % mu_basis: Basis column value for the solution vector

    [d, m] = size(matrix);
    %%% ! Step 1: Choose the basis column and initialize
    active_cols = find(any(matrix, 1)); % Find all columns with nonzero entries
    basis_col = min(active_cols); % Choose the smallest column as the basis
    gap_vector = -inf(1, m); % Initialize the gap vector with -inf
    gap_vector(basis_col) = 0; % The gap for the basis column is 0
    
    %%% ! Do search based on the basis column
    % Step 2: Use Breadth-First Search (BFS) to compute gaps
    queue = basis_col; % Initialize the queue with the basis column
    while ~isempty(queue)
        current_col = queue(1); % Take the first column from the queue
        queue(1) = []; % Remove it from the queue
        
        % Find all rows where the current column is active (nonzero)
        active_rows = find(matrix(:, current_col) ~= 0);
        %%% ! Important step: active_rows are those with nonzero entries
        % Iterate over these rows
        for i = 1:length(active_rows)
            row_idx = active_rows(i);
            
            % Find all columns with nonzero entries in this row
            other_active_cols = find(matrix(row_idx, :) ~= 0);
            
            % Iterate over these columns
            for j = 1:length(other_active_cols)
                col_idx = other_active_cols(j);
                
                %%% ! If this column hasn't been used for search yet
                if gap_vector(col_idx) == -inf
                    % Compute the gap: current column's gap + (current row's column value - current row's basis column value)
                    gap_vector(col_idx) = gap_vector(current_col) + ...
                                          (matrix(row_idx, col_idx) - matrix(row_idx, current_col));
                    
                    % Add this column to the queue for further processing
                    queue = [queue, col_idx];
                end
            end
        end
    end
    
    % Step 3: Compute mu_basis
    mu_basis = log(sum(B_subset)) - log(sum(exp(gap_vector(gap_vector ~= -inf))));
end