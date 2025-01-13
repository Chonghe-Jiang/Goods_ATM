function classes = quasi_equi(matrix)
    % Input: matrix - a binary matrix where matrix(i,j) = 1 indicates a connection between row i and column j
    % Output: classes - a cell array where each cell contains the indices of rows that belong to the same class

    [numRows, numCols] = size(matrix); % Get the dimensions of the matrix
    visited = false(1, numRows); % Track visited rows
    classes = {}; % Store the classes

    for i = 1:numRows
        if ~visited(i)
            % Start a new class with the current row
            class = [];
            stack = i;

            % Perform Depth-First Search (DFS) to find all connected rows
            while ~isempty(stack)
                row = stack(end);
                stack(end) = [];
                
                if ~visited(row)
                    visited(row) = true;
                    class = [class, row];
                    
                    % Find all columns connected to the current row
                    connectedCols = find(matrix(row, :) == 1);
                    
                    % Find all rows connected to these columns
                    for col = connectedCols
                        connectedRows = find(matrix(:, col) == 1)';
                        stack = [stack, setdiff(connectedRows, class)];
                    end
                end
            end

            % Add the class to the list of classes
            classes{end+1} = class;
        end
    end
end