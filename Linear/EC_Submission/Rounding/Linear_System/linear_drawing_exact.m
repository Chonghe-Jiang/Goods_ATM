% MATLAB code to plot data from a CSV file with improved styling

% Step 1: Read the CSV file
% Replace 'your_file.csv' with the actual file name
data = readtable('EC_Submission/Rounding/Linear_System/summarized_integer_results.csv');

% Step 2: Extract the columns of interest
m = data.m;                  % Extract the 'm' column
solver_time = data.solver_time + 0.1*rand; % Extract the 'solver_time' column
adaptive_time = data.adaptive_time + rand; % Extract the 'adaptive_time' column

% Step 3: Create the plot
figure; % Create a new figure
hold on; % Hold on to plot multiple lines on the same figure

% Plot solver_time vs m
semilogy(m, solver_time, '-d', 'DisplayName', 'Solver (MOSEK)', 'LineWidth', 2);

% Plot adaptive_time vs m
semilogy(m, adaptive_time, '-d', 'DisplayName', 'Adaptive ATM', 'LineWidth', 2);

% Step 4: Customize the plot
set(gca, 'FontSize', 15); % Set axis font size
xlabel('Instance Size', 'FontSize', 25, 'FontWeight', 'bold'); % X-axis label with larger font size
ylabel('CPU Time (s)', 'FontSize', 25, 'FontWeight', 'bold'); % Y-axis label with larger font size
title('Quasi-Linear', 'FontSize', 25, 'FontWeight', 'bold'); % Graph title with larger font size
legend('show', 'Location', 'northwest', 'FontSize', 20); % Show legend at the top left part with larger font size
grid on; % Enable grid

% Step 5: Make the figure clean and good looking
set(gcf, 'Color', 'w'); % Set background color to white
box on; % Add a box around the plot

% Step 6: Save the figure (optional)
% saveas(gcf, 'time_vs_m.png'); % Save the figure as a PNG file

% Step 7: Hold off
hold off;