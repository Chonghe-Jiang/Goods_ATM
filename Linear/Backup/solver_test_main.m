% 示例用法
n = 3;  % 约束数量
m = 4;  % 变量p的数量
v = rand(n, m);  % 随机生成v_ij
B = rand(n, 1);  % 随机生成B_i

res =  solver_test(v, B, n, m);

disp('Optimal p:');
disp(p_opt);
disp('Optimal beta:');
disp(beta_opt);