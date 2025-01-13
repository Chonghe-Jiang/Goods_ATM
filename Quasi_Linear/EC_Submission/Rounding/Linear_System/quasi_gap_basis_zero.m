function [gap_vector, mu_basis] = quasi_gap_basis_zero(matrix, B_subset)
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
    
    %%% ! Step 3: Compute mu_basis in a new way -> it is -a_0 now
    mu_basis = -gap_vector(1);
end