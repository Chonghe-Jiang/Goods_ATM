clc;
clear;

%%% Define the CSV filename
csv_filename = 'solver_adaptive_gap_data_repeat.csv';

%%% Read the CSV file into a table
data_table = readtable(csv_filename);

%%% Get unique classes based on the first two columns (n and m)
unique_classes = unique(data_table(:, 1:2), 'rows');

%%% Initialize a table to store the summarized results
summary_table = table('Size', [size(unique_classes, 1), width(data_table)], ...
    'VariableTypes', varfun(@class, data_table, 'OutputFormat', 'cell'), ...
    'VariableNames', data_table.Properties.VariableNames);

%%% Loop through each unique class
for i = 1:size(unique_classes, 1)
    % Extract the current class (n and m)
    current_n = unique_classes.n(i);
    current_m = unique_classes.m(i);
    
    % Find all rows in the data table that belong to this class
    class_rows = data_table(data_table.n == current_n & data_table.m == current_m, :);
    
    % Calculate the average of each column for this class
    avg_solver_time = mean(class_rows.solver_time);
    avg_adaptive_time = mean(class_rows.adaptive_time);
    avg_gap_optimal = mean(class_rows.gap_optimal);
    avg_Exact_solver_gap = mean(class_rows.Exact_solver_gap);
    avg_total_iter_adaptive = mean(class_rows.total_iter_adaptive);
    
    % Store the summarized results in the summary table
    summary_table.n(i) = current_n;
    summary_table.m(i) = current_m;
    summary_table.solver_time(i) = avg_solver_time;
    summary_table.adaptive_time(i) = avg_adaptive_time;
    summary_table.gap_optimal(i) = avg_gap_optimal;
    summary_table.Exact_solver_gap(i) = avg_Exact_solver_gap;
    summary_table.total_iter_adaptive(i) = avg_total_iter_adaptive;
end

%%% Print the summarized results
disp('Summarized Results:');
disp(summary_table);

%%% Optionally, save the summarized results to a new CSV file
%%% ! Summarize results to a specified CSV file - integer/exponential/rand
summary_csv_filename = 'summarized_integer_results.csv';
writetable(summary_table, summary_csv_filename);
disp(['Summarized results saved to ', summary_csv_filename]);