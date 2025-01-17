function classes = quasi_equi_optimized(matrix)
    % Input: matrix - a binary matrix where matrix(i,j) = 1 indicates a connection between row i and column j
    % Output: classes - a cell array where each cell contains the indices of rows that belong to the same class
    [numRows, numCols] = size(matrix); % Get the dimensions of the matrix
    visited = false(1, numRows); % Track visited rows
    classes = {}; % Store the classes

    % Precompute row-to-column and column-to-row connections
    rowToCols = cell(numRows, 1); % Each cell contains the columns connected to the row
    colToRows = cell(numCols, 1); % Each cell contains the rows connected to the column

    for i = 1:numRows
        rowToCols{i} = find(matrix(i, :) == 1); % Columns connected to row i
    end

    for j = 1:numCols
        colToRows{j} = find(matrix(:, j) == 1)'; % Rows connected to column j
    end

    % Preallocate stack and stack_top
    stack = zeros(1, numRows); % Preallocate stack with maximum possible size
    stack_top = 0; % Track the top of the stack

    for i = 1:numRows
        if ~visited(i)
            % Start a new class with the current row
            class = [];
            stack_top = stack_top + 1; % Push current row to stack
            stack(stack_top) = i;

            % Perform Depth-First Search (DFS) to find all connected rows
            while stack_top > 0
                % Pop the top element from the stack
                row = stack(stack_top);
                stack_top = stack_top - 1;

                if ~visited(row)
                    visited(row) = true;
                    class = [class, row];

                    % Find all columns connected to the current row
                    connectedCols = rowToCols{row};

                    % Find all rows connected to these columns
                    for col = connectedCols
                        connectedRows = colToRows{col};
                        for r = connectedRows
                            if ~visited(r)
                                stack_top = stack_top + 1; % Push unvisited row to stack
                                stack(stack_top) = r;
                            end
                        end
                    end
                end
            end

            % Add the class to the list of classes
            classes{end+1} = class;
        end
    end
end