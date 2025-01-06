n = 3;  % Number of rows
m = 3;  % Number of columns
A = randi([0, 1], n, m);

% Display the original matrix A
disp('Original Matrix A:');
disp(A);

% Call the function to find submatrices
result = findSubmatrices(A);

% Display the result
disp('Result:');
disp(result);

%%% Todo: Complete this function
%%% Todo: Complete the rounding procedure
