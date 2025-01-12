function [solution_adaptive, total_time_adaptive, total_iter_adaptive, obj_adaptive, dis_adaptive] = linear_dual_adaptive(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, adaptive_plot_flag, plot_flag_smooth,  p_opt_solver, fval_solver, adaptive, phase_num)

%%% Todo - stopping/switching
time_adaptive = [];
iter_adaptive = [];
obj_adaptive = [];
dis_adaptive = [];
total_time_adaptive = 0;  % Initialize total time
total_iter_adaptive = 0;
% Loop over phases
for phase = 1:phase_num
    % Call the linear_dual_agd function
    [solution_phase, time_phase, iter_phase, obj_phase, dis_phase, convergence] = linear_dual_agd(v, B, mu_0, max_iter, L, sigma, epsilon, mu_lower, mu_upper, delta, plot_flag, plot_flag_smooth, p_opt_solver, fval_solver, adaptive);
    
    % Update variables
    time_adaptive = [time_adaptive; time_phase];
    iter_adaptive = [iter_adaptive; iter_phase];
    obj_adaptive = [obj_adaptive; obj_phase];
    dis_adaptive = [dis_adaptive; dis_phase];
    
    % Accumulate total time
    total_time_adaptive = total_time_adaptive + time_phase;
    total_iter_adaptive = total_iter_adaptive + iter_phase;
    
    % Print phase information
    fprintf('Phase %d: Iterations = %d, Initial Objective = %f, Final Objective = %f\n', phase, iter_phase, obj_phase(1), obj_phase(end));
    %%% Todo - note that the stopping criteria (i.e. convergence) and the switching mechanism have not been available
    % Check for convergence
    if convergence
        break;
    end
    
    % Update delta and L
    %%% Todo - note that the switching way has not been determined
    if phase >=1 && phase <=3
        delta = delta / 3; 
        L = exp(max(mu_upper)) + sum(B) / delta;
        mu_0 = solution_phase;
    elseif phase <= 5
        delta = delta / 3; % Todo: previous 2 - for real data
        L = exp(max(mu_upper)) + sum(B) / delta;
        mu_0 = solution_phase;
    else
        delta = delta / 1.5;
        L = exp(max(mu_upper)) + sum(B) / delta;
        mu_0 = solution_phase;
    end
    % Return the final solution and total time
    solution_adaptive = mu_0;
end
% Plot the results if plot_flag is true
if adaptive_plot_flag
    figure;
    iterations = 1:length(obj_adaptive);
    
    subplot(2, 1, 1);
    plot(iterations, abs(obj_adaptive), 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Function Value Gap');
    title('Adaptive SAG - Function Value Convergence');
    grid on;
    
    subplot(2, 1, 2);
    plot(iterations, dis_adaptive, 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Iteration Distance');
    title('Adaptive SAG - Iteration Convergence');
    grid on;
end
end