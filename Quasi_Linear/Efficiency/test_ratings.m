clc
clear
% Assume your matrix is A
A = readmatrix('C:\Users\s1155203585\Desktop\EG_EXP\Linear\v0_funding_supp\Efficiency\Ratings_kroer.csv');

% Check if all entries are greater than or equal to 0
all_non_negative = all(A(:) >= 0);

% Display the result
if all_non_negative
    disp('All entries are greater than or equal to 0.');
else
    disp('Not all entries are greater than or equal to 0.');
end

% Check if the matrix contains any complex elements
if any(~isreal(A), 'all')
    disp('The matrix contains complex elements.');
else
    disp('The matrix does not contain complex elements.');
end

% Check if the matrix contains any NaN elements
if any(isnan(A), 'all')
    disp('The matrix contains NaN elements.');
else
    disp('The matrix does not contain NaN elements.');
end

% Check if the matrix contains any Inf elements
if any(isinf(A), 'all')
    disp('The matrix contains Inf elements.');
else
    disp('The matrix does not contain Inf elements.');
end

% Statistical survey of the matrix
min_value = min(A(:));
max_value = max(A(:));
mean_value = mean(A(:));
std_value = std(A(:));
matrix_size = size(A);

% Display the statistical results
fprintf('Matrix size: %dx%d\n', matrix_size(1), matrix_size(2));
fprintf('Minimum value: %.4f\n', min_value);
fprintf('Maximum value: %.4f\n', max_value);
fprintf('Mean value: %.4f\n', mean_value);
fprintf('Standard deviation: %.4f\n', std_value);