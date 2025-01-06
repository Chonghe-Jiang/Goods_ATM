function [p_opt, beta_opt] = solve_nonlinear_problem(v, B, n, m)
    % 检查输入维度
    assert(n > 0 && m > 0, 'n and m must be positive');
    assert(size(v, 1) == n && size(v, 2) == m, 'v must be n x m');
    assert(length(B) == n, 'B must have length n');

    % 初始化MOSEK问题结构
    prob = [];
    
    % 设置变量
    numvar = m + n;  % p (m个), beta (n个)
    prob.c = [ones(m,1); -B];  % 目标函数中包含 p 和 -B*log(beta)
    
    % 设置约束 p_j >= v_ij * beta_i
    A = sparse(n*m, numvar);
    for i = 1:n
        for j = 1:m
            A((i-1)*m + j, j) = 1;
            A((i-1)*m + j, m+i) = -v(i,j);
        end
    end
    prob.a = A;
    prob.blc = zeros(n*m, 1);
    prob.buc = inf(n*m, 1);
    
    % 设置变量界限
    prob.blx = [-inf(m,1); zeros(n,1)];  % beta 非负，p 无界
    prob.bux = inf(numvar, 1);
    
    % 设置指数锥约束 exp(-t_i/B_i) <= beta_i
    prob.cones = cell(n, 1);
    for i = 1:n
        prob.cones{i}.type = 'PEXP';
        prob.cones{i}.sub = [m+i, 1, 1];  % 这里需要确保维度正确
    end
    
    % 求解问题
    [~, res] = mosekopt('minimize echo(0)', prob);
    
    % 检查求解结果
    if isfield(res, 'rcode') && res.rcode ~= 0
        error('MOSEK error: %s (code %d)', res.rmsg, res.rcode);
    end
    
    % 提取结果
    if isfield(res, 'sol') && isfield(res.sol, 'itr')
        p_opt = res.sol.itr.xx(1:m);
        beta_opt = res.sol.itr.xx(m+1:m+n);
    else
        error('Unexpected MOSEK output structure');
    end
end