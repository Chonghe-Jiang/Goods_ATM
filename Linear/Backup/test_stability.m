clc;
clear;
n = 500;
m = 500;
delta = 1e-3;
v = rand(n,m);
B = rand(n,1);
mu = rand(1,m);
f_smooth = sum(exp(mu)) + delta * sum(B .* log(sum(exp((log(v) - repmat(mu,n,1)) / delta), 2)))

%%% ! For comparison
f_original = sum(exp(mu)) + sum(B .* max(log(v)-mu, [], 2))

f_original_upper = f_original + delta*log(m)*sum(B)

%%% ! Check the new formulation
max_log_v_mu = max(log(v) - repmat(mu, n, 1), [], 2);
% Rescale the values by subtracting the minimum
rescaled_log_v_mu = (log(v) - repmat(mu, n, 1) - max_log_v_mu) / delta;

% Compute the log-sum-exp term using the rescaled values
log_sum_exp_term = log(sum(exp(rescaled_log_v_mu), 2));

% Compute the final expression
f_smooth_check = sum(exp(mu)) + delta * sum(B .* ((max_log_v_mu/delta) + log_sum_exp_term));

