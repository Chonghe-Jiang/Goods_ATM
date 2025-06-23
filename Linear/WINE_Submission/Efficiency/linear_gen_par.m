function [x0, p0, mu_0, max_iter, step_size, eta, epsilon, L, sigma, mu_lower, mu_upper, delta] = linear_gen_par(v, B)
    % Generate parameters based on input v and B
    
    % Dimensions
    [n, m] = size(v);
    
    % Maximum number of iterations
    max_iter = 4000;
    
    % Step size for the subgradient method, candidates include: 2/4/6 * e-5 (in kroer paper)
    % Step size for the mirror descent method
    step_size = 1e-4; 
    eta = 0.02;
    
    % Tolerance for convergence
    epsilon = 1e-2;
    
    %%% Todo: Just warning here
    % Lower bound for mu - remember the log
    mu_lower = log(max(v .* B ./ sum(abs(v))));
    
    % Upper bound for mu - remember the log
    mu_upper = log(norm(B, 1) * ones(1, m));
    
    % Initializatin
    P = @(mu) max(mu_lower, min(mu, mu_upper));
    % mu_0 = P(log(1/m * ones(1,m)));
    % mu_0 = (mu_lower + mu_upper) / 2;
    % p0 = exp(mu_0);
    % x0 = rand(n,m);
    % Generate x0 in the middle of the area
    x0 = 1/m * ones(n, m).*B;
    % temp
    p0 = sum(x0,1);
    mu_0 = log(p0);
    % Final
    mu_0 = P(mu_0);
    p0 = exp(mu_0);
    
    % Delta parameter for the smoothing parameter - normally we do not have
    delta = 0.01;
    
    % Strong convexity parameter
    sigma = min(exp(mu_lower));
    
    % Lipschitz constant - its is a conservative estimate of the lipschitz constant
    L = (exp(max(mu_upper)) + sum(B) / delta); % without 2 version
    
    % % Print parameters to command line
    % fprintf('Generated Parameters:\n');
    % fprintf('p0:            [%s]\n', sprintf(' %f', p0));
    % fprintf('mu_0:          [%s]\n', sprintf(' %f', mu_0));
    % fprintf('max_iter:      %d\n', max_iter);
    % fprintf('step_size_subgradient:     %f\n', step_size);
    % fprintf('epsilon:       %e\n', epsilon);
    % fprintf('Lipschitz:             %f\n', L);
    % fprintf('SC:         %f\n', sigma);
    % fprintf('Smoothing term:         %f\n', delta);
end