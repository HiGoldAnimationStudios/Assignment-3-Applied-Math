function plot_multi_h_ref()
    tspan=[0,10];
    X0=1;
    h_ref=0.1;
    [t_list_fe,X_list_fe,~, ~] = fixed_step_integration(@rate_func01,@forward_euler_step,tspan,X0,h_ref);
    [t_list_be,X_list_be,~, ~] = fixed_step_integration(@rate_func01,@backward_euler_step,tspan,X0,h_ref);
    [t_list_me,X_list_me,~, ~] = fixed_step_integration(@rate_func01,@explicit_midpoint_step,tspan,X0,h_ref);
    [t_list_mi,X_list_mi,~, ~] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref);

    h_ref=0.2;
    [t_list_fe2,X_list_fe2,~, ~] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_be2,X_list_be2,~, ~] = fixed_step_integration(@rate_func01,@backward_euler_step,tspan,X0,h_ref);
    [t_list_me2,X_list_me2,~, ~] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_mi2,X_list_mi2,~, ~] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref);

    h_ref=0.3;
    [t_list_fe3,X_list_fe3,~, ~] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_be3,X_list_be3,~, ~] = fixed_step_integration(@rate_func01,@backward_euler_step,tspan,X0,h_ref);
    [t_list_me3,X_list_me3,~, ~] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_mi3,X_list_mi3,~, ~] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref);

    h_ref=0.4;
    [t_list_fe4,X_list_fe4,~, ~] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_be4,X_list_be4,~, ~] = fixed_step_integration(@rate_func01,@backward_euler_step,tspan,X0,h_ref);
    [t_list_me4,X_list_me4,~, ~] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_mi4,X_list_mi4,~, ~] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref);
    

    
    t = linspace(0, 10, 400);
    % Forward Euler
    figure(); clf; hold on; grid on;
    plot(t, solution01(t), 'k', 'LineWidth', 1.5, 'DisplayName','Exact');
    plot(t_list_fe,  X_list_fe,  'o-', 'Color',[1 0 0], 'MarkerSize',5, 'LineWidth', 0.5, 'DisplayName','Forward Euler  h=0.1');
    plot(t_list_fe2, X_list_fe2, 'o-', 'Color',[0 0 1], 'MarkerSize',5, 'LineWidth', 0.5, 'DisplayName','Forward Euler  h=0.2');
    plot(t_list_fe3, X_list_fe3, 'o-', 'Color',[0 0.5 0], 'MarkerSize',5, 'LineWidth', 0.5, 'DisplayName','Forward Euler  h=0.3');
    plot(t_list_fe4, X_list_fe4, 'o-', 'Color',[0.85 0.65 0], 'MarkerSize',5, 'LineWidth', 0.5, 'DisplayName','Forward Euler  h=0.4');
    xlabel('t'); ylabel('X(t)');
    title('Forward Euler vs. Exact Solution');
    legend('Location','best'); hold off;
    
    % Explicit Midpoint
    figure(); clf; hold on; grid on;
    plot(t, solution01(t), 'k-', 'LineWidth',1.5, 'DisplayName','Exact');
    plot(t_list_me,  X_list_me,  'o-', 'Color',[1 0 0],   'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Explicit Midpoint  h=0.1');
    plot(t_list_me2, X_list_me2, 'o-', 'Color',[0 0 1],   'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Explicit Midpoint  h=0.2');
    plot(t_list_me3, X_list_me3, 'o-', 'Color',[0 0.5 0], 'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Explicit Midpoint  h=0.3');
    plot(t_list_me4, X_list_me4, 'o-', 'Color',[0.85 0.65 0], 'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Explicit Midpoint  h=0.4');
    xlabel('t'); ylabel('X(t)'); title('Explicit Midpoint vs. Exact Solution');
    legend('Location','best'); hold off;
    
    % Backward Euler
    figure(); clf; hold on; grid on;
    plot(t, solution01(t), 'k-', 'LineWidth',1.5, 'DisplayName','Exact');
    plot(t_list_be,  X_list_be,  'o-', 'Color',[1 0 0],   'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Backward Euler  h=0.1');
    plot(t_list_be2, X_list_be2, 'o-', 'Color',[0 0 1],   'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Backward Euler  h=0.2');
    plot(t_list_be3, X_list_be3, 'o-', 'Color',[0 0.5 0], 'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Backward Euler  h=0.3');
    plot(t_list_be4, X_list_be4, 'o-', 'Color',[0.85 0.65 0], 'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Backward Euler  h=0.4');
    xlabel('t'); ylabel('X(t)'); title('Backward Euler vs. Exact Solution');
    legend('Location','best'); hold off;
    
    % Implicit Midpoint 
    figure(); clf; hold on; grid on;
    plot(t, solution01(t), 'k-', 'LineWidth',1.5, 'DisplayName','Exact');
    plot(t_list_mi,  X_list_mi,  'o-', 'Color',[1 0 0],   'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Implicit Midpoint  h=0.1');
    plot(t_list_mi2, X_list_mi2, 'o-', 'Color',[0 0 1],   'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Implicit Midpoint  h=0.2');
    plot(t_list_mi3, X_list_mi3, 'o-', 'Color',[0 0.5 0], 'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Implicit Midpoint  h=0.3');
    plot(t_list_mi4, X_list_mi4, 'o-', 'Color',[0.85 0.65 0], 'MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Implicit Midpoint  h=0.4');
    xlabel('t'); ylabel('X(t)'); title('Implicit Midpoint vs. Exact Solution');
    legend('Location','best'); hold off;
end

function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end
