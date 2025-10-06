function day10()
    tspan=[0,10];
    X0=1;
    h_ref=0.1;
    [t_list_e,X_list_e,h_avg_e, num_evals_e] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_m,X_list_m,h_avg_m, num_evals_m] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    h_ref=0.2;
    [t_list_e2,X_list_e2,h_avg_e2, num_evals_e2] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    [t_list_m2,X_list_m2,h_avg_m2, num_evals_m2] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);

    t = linspace(0, 10, 400);
    plot(t, solution01(t),"k"), grid on
    hold on
    plot(t_list_e, X_list_e,"r--")
    plot(t_list_m, X_list_m,"b--")
    plot(t_list_e2, X_list_e2,"r")
    plot(t_list_m2, X_list_m2,"b")
    
    t_ref=0.3676767676767;
    local_trunc_error(@rate_func01, @solution01, t_ref, X0, h_ref);
    
end

function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end