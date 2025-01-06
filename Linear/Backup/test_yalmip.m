% Clear workspace and command window
clear;
clc;

% Define problem parameters
n = 5; % Range of i
m = 3; % Range of j
B = rand(n, 1); % Random values for B_i
v = rand(n, m); % Random values for v_{ij}

% Define variables
p = sdpvar(m, 1); % p_j variables
beta = sdpvar(n, 1); % beta_i variables

% Define objective function
objective = sum(p) - sum(B .* log(beta));

% Define constraints
constraints = [];
for i = 1:n
    for j = 1:m
        constraints = [constraints, p(j) >= v(i, j) * beta(i)];
    end
end

% Set optimization problem
options = sdpsettings('solver', 'mosek');

% Solve optimization problem
sol = optimize(constraints, objective, options);

% Check solution status
if sol.problem == 0
    % Extract solutions
    p_opt = value(p);
    beta_opt = value(beta);
    disp('Optimization successful');
    disp('Optimal values for p:');
    disp(p_opt);
    disp('Optimal values for beta:');
    disp(beta_opt);
else
    disp('Optimization failed');
    disp('Error message:');
    disp(sol.info);
end