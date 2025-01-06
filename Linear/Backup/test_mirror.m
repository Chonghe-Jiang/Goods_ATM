% test_mirror.m
% Contents: test the correctness of our mirror descent algorithm
% Using the distance in every iteration between optimal and the mirror iterates

% Generate random v matrix and B vector
n = 10; % number of rows
m = 10; % number of columns

% generate v and B
% Generate budgets from uniform distribution and normalize
B = rand(n, 1); % Draw budgets from uniform distribution
B = B / sum(B); % Normalize to sum to 1

%% solvers - optimal 

% Generate valuations from exponential distribution and normalize
v = exprnd(1, n, m ); % Draw valuations from exponential distribution
v = v ./ sum(v, 2); % Normalize each row to sum to 1
plot_flag = true;

% Parameters for the algorithms
[x0, p0, mu_0, max_iter, step_size, eta, epsilon, L, sigma, mu_lower, mu_upper, delta] = linear_gen_par(v, B);

% Test linear_dual_solver
fprintf('Testing linear_dual_solver...\n');
[p_opt_solver, fval_solver] = linear_dual_solver(v, B, p0);
fprintf('Optimal function value from solver is: %.4f\n', fval_solver);

% Test the correctness of mirror descent
x_mirror = linear_primal_md(v, B, x0, eta, epsilon, max_iter, p_opt_solver, plot_flag);
p_mirror = sum(x_mirror, 1);
obj_mirror = sum(p_mirror) - sum(B .* log(min(p_mirror ./ v, [], 2)));
fprintf('Optimal function value from MD is: %.4f\n', obj_mirror);

