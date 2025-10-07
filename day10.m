function day10()
    tspan=[0,10];
    X0=1;
    h_ref=0.1;
    [t_list_fe,X_list_fe,h_avg_fe, num_evals_fe] = fixed_step_integration(@rate_func01,@forward_euler_step,tspan,X0,h_ref);
    [t_list_be,X_list_be,h_avg_be, num_evals_be] = fixed_step_integration(@rate_func01,@backward_euler_step,tspan,X0,h_ref);
    [t_list_me,X_list_me,h_avg_me, num_evals_me] = fixed_step_integration(@rate_func01,@explicit_midpoint_step,tspan,X0,h_ref);
    [t_list_mi,X_list_mi,h_avg_mi, num_evals_mi] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref);

    h_ref=0.2;
    [t_list_fe2,X_list_fe2,h_avg_fe2, num_evals_fe2] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_be2,X_list_be2,h_avg_be2, num_evals_be2] = fixed_step_integration(@rate_func01,@backward_euler_step,tspan,X0,h_ref);
    [t_list_me2,X_list_me2,h_avg_me2, num_evals_me2] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_mi2,X_list_mi2,h_avg_mi2, num_evals_mi2] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref);

    t = linspace(0, 10, 400);
    figure(1);
    plot(t, solution01(t),"k"), grid on
    hold on
    plot(t_list_fe, X_list_fe,"r--")
    plot(t_list_fe2, X_list_fe2,"r")
    hold off
    
    figure(2);
    plot(t, solution01(t),"k"), grid on
    hold on
    plot(t_list_me, X_list_me,"b--")
    plot(t_list_me2, X_list_me2,"b")
    hold off

    figure(3);
    plot(t, solution01(t),"k"), grid on
    hold on
    plot(t_list_be, X_list_be,"g--")
    plot(t_list_be2, X_list_be2,"g")
    hold off

    figure(4);
    plot(t, solution01(t),"k"), grid on
    hold on
    plot(t_list_mi, X_list_mi,"y--")
    plot(t_list_mi2, X_list_mi2,"y")
    hold off

    t_ref = 0.3676767676767;

    local_trunc_error(@rate_func01, @solution01, t_ref);
    
end

function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end