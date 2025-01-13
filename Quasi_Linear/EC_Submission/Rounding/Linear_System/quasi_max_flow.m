function [is_feasible, optimal_value, sum_exp_mu] = quasi_max_flow(mu, B, v, R)
    % quasi_max_flow_extended: Solves the max-flow problem with an extra dimension
    % Inputs:
    %   mu: 1 x (m+1) vector
    %   B: n x 1 vector
    %   R: n x (m+1) binary matrix
    %   v: n x (m+1) matrix
    % Outputs:
    %   is_feasible: True if the optimal value equals sum(exp(mu)), False otherwise
    %   optimal_value: Optimal value of the max-flow problem
    %   sum_exp_mu: Sum of exp(mu)
    % ! Already with updated input
    % Dimensions
    m_plus_1 = length(mu); % m+1
    n = length(B);
    % Step 1: Check if linear_activation(mu, v) is equal to R
    if ~isequal(linear_activation(mu, v), R)
        is_feasible = false;
        optimal_value = inf;
        sum_exp_mu = inf;
        return;
    end
    
    V = {'s', 't'};
    for j = 0:m_plus_1-1
        V{end+1} = sprintf('g%d', j);
    end
    for i = 1:n
        V{end+1} = sprintf('b%d', i);
    end

    % Define the edge set E and capacities
    E = [];
    C = [];

    % Edges from source 's' to g_j with capacity exp(mu_j)
    for j = 0:m_plus_1-1
        E = [E; {'s', sprintf('g%d', j)}];
        C = [C; exp(mu(j+1))]; % Adjust index for MATLAB (1-based)
    end

    % Edges from g_j to b_i with capacity +inf if R(i,j+1) == 1
    for i = 1:n
        for j = 0:m_plus_1-1
            if R(i,j+1) == 1 % Adjust index for MATLAB (1-based)
                E = [E; {sprintf('g%d', j), sprintf('b%d', i)}];
                C = [C; inf];
            end
        end
    end

    % Edges from b_i to sink 't' with capacity B_i
    for i = 1:n
        E = [E; {sprintf('b%d', i), 't'}];
        C = [C; B(i)];
    end

    % Create the graph
    G = digraph(E(:,1), E(:,2), C);

    % Solve the max-flow problem
    mf = maxflow(G, 's', 't');

    optimal_value = mf;
    sum_exp_mu = sum(exp(mu));
    
    % Check feasibility with numerical precision
    is_feasible = abs(optimal_value - sum_exp_mu) < 1e-2;
end

function binary_matrix = linear_activation(mu, v)
    % linear_activation - Given a 1 x (m+1) vector mu and an n x (m+1) matrix v, output a binary matrix
    % where each element is 1 if it is the maximum of its row, and 0 otherwise.
    %
    % Input:
    % mu - 1 x (m+1) vector
    % v - n x (m+1) matrix
    %
    % Output:
    % binary_matrix - n x (m+1) binary matrix

    % Calculate the adjusted values
    adjusted_values = log(v) - mu;

    % Find the maximum element in each row
    max_values = max(adjusted_values, [], 2);

    % Create a binary matrix where each element is 1 if it is the maximum of its row
    % Use a numerical precision threshold of 1e-8
    binary_matrix = abs(adjusted_values - max_values) < 1e-6; %%% Todo: Need to revise this
end