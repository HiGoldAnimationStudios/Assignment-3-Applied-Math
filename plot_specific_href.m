function plot_specific_href(h)
    % Setup
    tspan = [0, 20];
    X0    = 1;
    
    % exact
    t_plot = linspace(tspan(1), tspan(2), 600);
    x_exact = solution01(t_plot);

    % Numerical solutions
    [t_fe, X_fe] = fixed_step_integration(@rate_func01, @forward_euler_step,     tspan, X0, h);
    [t_em, X_em] = fixed_step_integration(@rate_func01, @explicit_midpoint_step, tspan, X0, h);
    [t_be, X_be] = fixed_step_integration(@rate_func01, @backward_euler_step,    tspan, X0, h);
    [t_im, X_im] = fixed_step_integration(@rate_func01, @implicit_midpoint_step, tspan, X0, h);

    figure('Color','w');
    sgtitle(sprintf('Numerical Solutions vs Exact Solution (h = %.2f)', h))

    % Forward Euler
    subplot(4,1,1); hold on; grid on;
    plot(t_plot, x_exact, 'k-', 'LineWidth',0.5, 'DisplayName','Exact');
    plot(t_fe,   X_fe,   'o-', 'LineWidth',0.5, 'MarkerSize',5, 'Color',[1 0 0], 'DisplayName','Forward Euler');
    xlabel('t'); ylabel('X(t)');
    title(sprintf('Forward Euler (h = %.2f)', h))
    legend('Location','best');

    % Explicit Midpoint
    subplot(4,1,2); hold on; grid on;
    plot(t_plot, x_exact, 'k-', 'LineWidth',0.5, 'DisplayName','Exact');
    plot(t_em,   X_em,   'o-', 'LineWidth',0.5, 'MarkerSize',5, 'Color',[0 0.5 0], 'DisplayName','Explicit Midpoint');
    xlabel('t'); ylabel('X(t)');
    title(sprintf('Explicit Midpoint (h = %.2f)', h))
    legend('Location','best');

    % Backward Euler
    subplot(4,1,3); hold on; grid on;
    plot(t_plot, x_exact, 'k-', 'LineWidth',0.5, 'DisplayName','Exact');
    plot(t_be,   X_be,   'o-', 'LineWidth',0.5, 'MarkerSize',5, 'Color',[0 0 1], 'DisplayName','Backward Euler');
    xlabel('t'); ylabel('X(t)');
    title(sprintf('Backward Euler (h = %.2f)', h))
    legend('Location','best');

    % Implicit Midpoint
    subplot(4,1,4); hold on; grid on;
    plot(t_plot, x_exact, 'k-', 'LineWidth',0.5, 'DisplayName','Exact');
    plot(t_im,   X_im,   'o-', 'LineWidth',0.5, 'MarkerSize',5, 'Color',[0.7 0 0.7], 'DisplayName','Implicit Midpoint');
    xlabel('t'); ylabel('X(t)');
    title(sprintf('Implicit Midpoint (h = %.2f)', h))
    legend('Location','best');
end

function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end
