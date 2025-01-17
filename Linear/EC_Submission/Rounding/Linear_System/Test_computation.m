matrix = [1,1,0;0,0,1;0,0,1];
% matrix = sparse(matrix);
% 测试优化后的函数
tic;
classes = linear_equi_optimized(matrix);
toc;

% 显示等价类的数量
disp(['Class number: ', num2str(length(classes))]);