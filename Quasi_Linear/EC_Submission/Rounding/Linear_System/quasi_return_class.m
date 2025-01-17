function [row_classes, column_classes, col_class_matrices] = quasi_return_class(A, v)
    % Input:
    % A - binary input matrix (rows represent nodes, columns represent features)
    % v - a matrix of the same size as A, used to compute log values
    % Output:
    % row_classes - a cell array where each cell contains the row indices of a class
    % column_classes - a cell array where each cell contains the column indices as a row vector
    % col_class_matrices - a cell array where each cell contains the derived matrix for a row class

    % Step 1: Find row classes (connected components) in A
    row_classes = quasi_equi_optimized(A);
    for i = 1:numel(row_classes)
        row_classes{i} = sort(row_classes{i});
    end

    % Initialize outputs
    col_class_matrices = cell(1, length(row_classes));
    column_classes = cell(1, length(row_classes));

    % Step 2: For each row class, derive the matrix and column class
    for i = 1:length(row_classes)
        same_rows = row_classes{i}; % Rows in the current row class
        num_rows = length(same_rows);
        num_cols = size(A, 2);

        % Initialize the matrix for this row class
        col_class = zeros(num_rows, num_cols); % Default to -inf (or any placeholder)

        % Initialize the column class for this row class
        column_class = [];

        % Fill the matrix and determine the column class
        for j = 1:num_rows
            row_idx = same_rows(j);
            active_cols = A(row_idx, :) == 1; % Columns where A(row_idx, :) == 1
            col_class(j, active_cols) = log(v(row_idx, active_cols)); % Apply log(v)
            column_class = union(column_class, find(active_cols)); % Add active columns to the column class
        end

        % Store the derived matrix and column class
        col_class_matrices{i} = col_class;
        column_classes{i} = column_class(:).'; % Ensure column class is a row vector
    end
end