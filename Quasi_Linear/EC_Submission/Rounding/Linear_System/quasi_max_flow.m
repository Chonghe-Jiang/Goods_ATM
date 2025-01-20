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
    % Step 1: Check if linear_activation(mu, v) is equal to R
    % ! This is for double check the machine accuracy
    sum(sum(abs(quasi_activation(mu, v)-R)))
    if ~isequal(quasi_activation(mu, v), R)
        is_feasible = false;
        optimal_value = inf;
        sum_exp_mu = inf;
        return;
    end
    mu = mu(2:end); % Remove the last element
    % 假设 mu, R, B, m, n 已经定义

    % 定义节点集合 V
    V = {'s', 't'};
    for j = 0:length(mu)  % 注意：j 从 0 开始
        V{end+1} = sprintf('g%d', j);
    end
    for i = 1:length(B)
        V{end+1} = sprintf('b%d', i);
    end

    % 定义边集 E 和容量 C
    E = [];
    C = [];

    % 从源节点 's' 到 g_j 的边，容量为 exp(mu_j)（j > 0）
    for j = 1:length(mu)  % 注意：j 从 1 开始，因为 j = 0 没有 mu_j
        E = [E; {'s', sprintf('g%d', j)}];
        C = [C; exp(mu(j))];
    end

    % j = 0 不需要连接到 mu_j，因此直接从 's' 到 g0 的边容量为无穷大
    E = [E; {'s', 'g0'}];
    C = [C; inf];

    % 从 g_j 到 b_i 的边，容量为 +inf 如果 R(i,j) == 1
    for i = 1:length(B)
        for j = 0:length(mu)  % 注意：j 从 0 开始
            if j == 0 || R(i,j) == 1  % j = 0 时总是连接，j > 0 时根据 R(i,j) 决定
                E = [E; {sprintf('g%d', j), sprintf('b%d', i)}];
                C = [C; inf];
            end
        end
    end

    % 从 b_i 到汇节点 't' 的边，容量为 B_i
    for i = 1:length(B)
        E = [E; {sprintf('b%d', i), 't'}];
        C = [C; B(i)];
    end

    % 创建有向图
    G = digraph(E(:,1), E(:,2), C);

    % 求解最大流问题
    mf = maxflow(G, 's', 't');

    optimal_value = mf;
    sum_exp_mu = sum(exp(mu));
    
    % Check feasibility with numerical precision
    is_feasible = abs(optimal_value - sum_exp_mu) < 1e-2;
end

function binary_matrix = quasi_activation(mu, v)
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
    binary_matrix = abs(adjusted_values - max_values) < 1e-4; %%% Todo: Need to revise this
end