function result = findSubmatrices(A)
    [n, m] = size(A);
    result = {};
    processed_rows = [];
    
    while any(A(:))
        % Find the next unprocessed row
        unprocessed_rows = setdiff(1:n, processed_rows);
        if isempty(unprocessed_rows)
            break;
        end
        i = unprocessed_rows(1);
        processed_rows = [processed_rows, i];
        
        P = i;  % Initialize P with the current row
        Q = [];
        
        % Find initial Q1
        Q1 = find(A(i, :) == 1);
        % if isempty(Q1)
        %     continue;
        % end
        Q = [Q, Q1];
        
        % Find initial P1
        P1 = find(any(A(:, Q1) == 1, 2));
        if numel(P1) == 1  % If P1 only contains the element i
            result{end+1} = {P, Q};
            A(P, Q) = 0;
            continue;  % Skip to the next iteration of the outer while loop
        end
        P = unique([P, P1']);  % Ensure P includes new elements
        
        % Iterate to find new P and Q
        while true
            new_Q = unique(find(any(A(P, :) == 1, 1)));
            if isempty(new_Q) || all(ismember(new_Q, Q))
                break;
            end
            Q = unique([Q, new_Q]);
            
            new_P = find(any(A(:, Q) == 1, 2));
            if isempty(new_P) || all(ismember(new_P, P))
                break;
            end
            P = unique([P, new_P']);
        end
        
        % Store the result
        result{end+1} = {P, Q};
        
        % Set the found submatrix to 0
        A(P, Q) = 0;
        
        % Update processed_rows with new_P
        processed_rows = unique([processed_rows, P]);
    end
end