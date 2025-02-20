function solution = quasi_class_generation(A, v, B)
    %%% * Given binary matrix - output the solution
    % A: Input matrix with n rows and m columns, entries are either 0 or 1
    % v: Constant matrix with the same size as A
    % B: n by 1 vector
    % solution: 1*m vector
    %%% ! Step 1: Do classification of different classes
    %%% Todo: the input A -> output classes; we use brand new function type to represent it 
    [row_classes, column_classes, col_class_matrices] = quasi_return_class(A, v);
    %%% ! Step 2: Do computation for every class - use the above information
    m = size(A, 2); %! Here the size is actually m+1
    solution = zeros(1, m); %! Change initialization vector since it is quasi-linear version
    %%% ! The change for quasi-linear utility
    for i = 1:length(row_classes)
        % Extract same_rows and col_class for the current class
        same_rows = row_classes{i}; % Rows in the current class
        col_class = col_class_matrices{i}; % Derived matrix for the current class
        B_class = B(same_rows, :);
        % Check if index 1 belongs to col_class{i}
        if ismember(1, column_classes{i}) % ! here we should use the column_classes instead of the matrix
            % Todo: Recheck this step
            % Use quasi_gap_basis_zero if index 1 is in col_class{i}
            [gap_vector, mu_basis] = quasi_gap_basis_zero(col_class, B_class);
        else
            % Otherwise, use quasi_gap_basis_normal
            [gap_vector, mu_basis] = quasi_gap_basis_normal(col_class, B_class);
        end
        
        % Proceed with the rest of the logic
        active_cols = find(gap_vector ~= -inf);
        solution(active_cols) = mu_basis + gap_vector(active_cols);
    end
    solution = solution(2:end); % ! Only choose the 2:end result
end
%%% ! This is a new approach for the gap calculation
