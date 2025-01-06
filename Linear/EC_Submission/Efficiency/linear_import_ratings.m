% Read the CSV file
data = readmatrix('Dataset\Ratings.csv');

% Display the size of the imported matrix
[rows, cols] = size(data);
fprintf('Imported matrix size: %d rows x %d columns\n', rows, cols);

% Display the first few rows and columns of the matrix
disp('First 5 rows and 5 columns of the imported matrix:');
disp(data(1:5, 1:5));