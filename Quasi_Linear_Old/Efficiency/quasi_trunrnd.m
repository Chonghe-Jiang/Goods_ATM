function A = quasi_trunrnd(n, m, condition_num)
    % Generate a random matrix of size n by m
    A = rand(n, m);
    
    % Find the largest and smallest elements in absolute value
    max_val = max(abs(A(:)));
    min_val = min(abs(A(:)));
    
    % Normalize the matrix so that the largest element is 1
    A = A / max_val;
    
    % Adjust the smallest element to achieve the desired condition number
    min_val = 1 / condition_num;
    
    % Scale the matrix to adjust the smallest element
    A = A * (1 - min_val) + min_val;
end