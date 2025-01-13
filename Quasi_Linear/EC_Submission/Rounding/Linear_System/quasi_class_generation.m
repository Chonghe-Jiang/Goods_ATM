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
    m = size(A, 2); % Number of columns in A
    solution = zeros(1, m); % Initialize solution vector
    %%% ! The change for quasi-linear utility
    for i = 1:length(row_classes)
        % Extract same_rows and col_class for the current class
        same_rows = row_classes{i}; % Rows in the current class
        col_class = col_class_matrices{i}; % Derived matrix for the current class
        B_class = B(same_rows, :);
        
        % Check if index 1 belongs to col_class{i}
        if ismember(1, col_class{i})
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
end
%%% ! This is a new approach for the gap calculation
