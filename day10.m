function day10()
    
    %Plots comparing the closed-form solution to equation 5 with the numerical approximation (for a few different time step sizes). The idea is that each plot should depict that numerical approximation of x(t) getting more and more accurate as the step size gets smaller and smaller.

    %plot_multi_h_ref()

    %Plots comparing the numerical solutions to eqn. 5 for a time step of href = .38 and a separate plot for href = .45 (as detailed in the explicit vs. implicit stability section).

    %plot_specific_href(0.38)
    
    %plot_specific_href(0.45)

    %local_trunc_error(@rate_func01, @solution01, 0.3);

    global_trunc_error(@rate_func01, @solution01, 0, 20);
    
end

function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end

