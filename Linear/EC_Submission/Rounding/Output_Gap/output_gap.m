%%% EC 2025 Submission
clc;
clear;
addpath('C:\Users\s1155203585\Dropbox\EG_EXP\Linear\EC_Submission\Efficiency');
% Define the matrix sizes to test
matrix_sizes = [20, 50, 100, 200];

% Define the distributions to test
distributions = {'normal', 'integer', 'exponential', 'lognormal'};

% Initialize a cell array to store the results
results = cell(length(distributions), length(matrix_sizes));

% Loop over each matrix size
for i = 1:length(matrix_sizes)
    n = matrix_sizes(i);
    m = matrix_sizes(i);
    
    % Generate budgets from uniform distribution and normalize
    B = rand(n, 1); % Draw budgets from uniform distribution
    % B = B / sum(B); % Normalize to sum to 1
    
    % Loop over each distribution
    for j = 1:length(distributions)
        dist = distributions{j};
        
        % Generate valuations based on the distribution
        switch dist
            case 'normal'
                v = rand(n, m); % Draw valuations from normal distribution
            case 'integer'
                v = randi([1, 10], n, m); % Draw valuations from integer distribution
            case 'exponential'
                v = exprnd(1, n, m); % Draw valuations from exponential distribution
            case 'lognormal'
                v = lognrnd(0, 1, n, m); % Draw valuations from log-normal distribution
        end
        
        % Normalize each row to sum to 1
        v = v ./ sum(v, 2);
        
        % Solve the optimization problem
        [p_opt, beta_opt, objective_value, solve_time] = linear_dual_solver(n, m, B, v);
        
        % Compute the gap
        mu = log(p_opt);
        radius = 0; % Todo: set radius=0 for complete output
        [gap, gap_array, matrix_backup] = linear_compute_gap(v, B, mu,radius);
        
        % Store the gap result
        results{j, i} = gap;
    end
end

% Create a table to store the results
result_table = cell2table(results, 'VariableNames', cellstr(num2str(matrix_sizes')), 'RowNames', distributions);

% Export the table to an Excel file
filename = 'gap_results.xlsx';
writetable(result_table, filename, 'WriteRowNames', true);

fprintf('Results have been saved to %s\n', filename);