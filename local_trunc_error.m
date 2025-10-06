function local_trunc_error(rate_func_in, X, t_ref, X0, h_ref)
    h_values = logspace(-5,1,100);
    
    e_fw_e = zeros(length(h_values));
    
    e_midpoint = zeros(size(h_values));
    sol = zeros(size(h_values));

    disp(length(h_values))
    disp(length(e_fw_e))

    for i = 1:length(h_values)
        h = h_values(i);

        sol = X(t_ref + h);

        xB_fw_e=forward_euler_step(rate_func_in,t_ref,X0,h_ref);
        e_fw_e(i) = abs(xB_fw_e - sol);
        xB_midpoint=explicit_midpoint_step(rate_func_in,t_ref,X0,h_ref);
        e_midpoint(i) = abs(xB_midpoint  - sol);

    end
    
    disp(length(h_values))
    disp(length(e_fw_e))

    loglog_fit(h, sol)
    hold on
    loglog_fit(h, e_fw_e)
    loglog_fit(h, e_midpoint)
    hold_off


end 