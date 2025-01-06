function solution = linear_equation_solution(row_classes, col_classes, B)
    % row_classes: Cell array of row index sets
    % col_classes: Cell array of column matrices
    % B: Constant vector of size n*1
    % solution: 1*m vector
    %%% Todo: check correctness
    %%% Todo: modify the path finding approach

    n = length(B);
    m = size(col_classes{1}, 2);
    solution = zeros(1, m);
    
    for i = 1:length(row_classes)
        row_class = row_classes{i};
        col_class = col_classes{i};
        
        % Find the active columns for this row class
        active_cols = find(any(col_class ~= 0, 1));
        
        % Initialize the gap vector
        gap = -inf(1, m);
        gap(active_cols) = 0;
        
        % Choose the first active column as j_fix
        j_fix = active_cols(1);
        
        % Compute the gap for each active column
        for j = 1:length(active_cols)
            col_idx = active_cols(j);
            if col_idx ~= j_fix
                % Find the path to connect col_idx and j_fix
                path = find_path(col_class, j_fix, col_idx);
                if ~isempty(path)
                    % Compute the gap
                    gap_value = sum(col_class(path(:, 1), path(:, 2)));
                    gap(col_idx) = gap_value;
                end
            end
        end
        
        % Compute fix_value
        sum_B = sum(B(row_class));
        sum_exp_gap = sum(exp(gap(active_cols)));
        fix_value = log(sum_B) - log(sum_exp_gap);
        
        % Update the solution vector
        solution(j_fix) = fix_value;
        for j = 1:length(active_cols)
            col_idx = active_cols(j);
            if col_idx ~= j_fix
                solution(col_idx) = fix_value + gap(col_idx);
            end
        end
    end
end

function path = find_path(col_class, j_fix, col_idx)
    % Find a path from col_idx to j_fix in the col_class matrix
    % This is a simple BFS-like search for demonstration purposes
    % In a real implementation, you might need a more sophisticated algorithm
    % to handle larger matrices and more complex paths
    
    [n, m] = size(col_class);
    visited = false(n, m);
    queue = [1, col_idx]; % Start from the first row and col_idx
    
    while ~isempty(queue)
        current = queue(1, :);
        queue(1, :) = [];
        
        if current(2) == j_fix
            path = current;
            return;
        end
        
        visited(current(1), current(2)) = true;
        
        % Explore neighbors
        for i = 1:n
            if ~visited(i, current(2)) && col_class(i, current(2)) ~= 0
                queue = [queue; i, current(2)];
            end
        end
        
        for j = 1:m
            if ~visited(current(1), j) && col_class(current(1), j) ~= 0
                queue = [queue; current(1), j];
            end
        end
    end
    
    path = []; % No path found
end