%Runs numerical integration using explicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in    during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_midpoint_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    %your code here
    %your code here
    t_list=[];
    X_list=[];
    XA=X0;
    t0=tspan(1);
    tf=tspan(2);

    N=ceil((tf-t0)/h_ref);
    h_avg=(tf-t0)/N;
    num_evals=0;

    for t = t0 : h_avg : tf
        t_list(end+1)=t;
        X_list(end+1)=XA;
        [XB,num_evals_i] = explicit_midpoint_step(rate_func_in,t,XA,h_avg);
        num_evals=num_evals+num_evals_i;
        XA=XB; 
    end
end