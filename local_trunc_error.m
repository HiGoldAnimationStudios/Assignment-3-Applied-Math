function local_trunc_error(rate_func_in, X, t_ref)
    h_values = logspace(-5,1,100);
    
    e_fw_e = zeros(size(h_values));
    e_midpoint = zeros(size(h_values));
    e_implicit_mid = zeros(size(h_values));
    e_backw_mid = zeros(size(h_values));

    analytical_value = zeros(size(h_values));




    %disp(length(h_values))
    %disp(length(e_fw_e))

    X0 = X(t_ref);

    for i = 1:length(h_values)
        h = h_values(i);

        analytical_value(i) = X(t_ref + h);

        xB_fw_e =forward_euler_step(rate_func_in,t_ref,X0,h);
        e_fw_e(i) = abs(xB_fw_e - analytical_value(i));

        xB_midpoint = explicit_midpoint_step(rate_func_in,t_ref,X0,h);
        e_midpoint(i) = abs(xB_midpoint  - analytical_value(i));

        xB_backw_e = backward_euler_step(rate_func_in,t_ref,X0,h);
        e_backw_e(i) = abs(xB_backw_e - analytical_value(i));

        xB_implic_midpoint = implicit_midpoint_step(rate_func_in,t_ref,X0,h);
        e_implic_midpoint(i) = abs(xB_implic_midpoint  - analytical_value(i));










    end
    
    disp(size(h_values))
    disp(size(analytical_value))
    disp(size(e_fw_e))
    disp(size(e_midpoint))

 
    [p_x_loc,k_x_loc] = loglog_fit(h_values, abs(analytical_value-X0))
    [p_fw_loc,k_fw_loc] =loglog_fit(h_values, e_fw_e)
    [p_ex_loc,k_ex_loc] =loglog_fit(h_values, e_midpoint)
    [p_im_loc,k_im_loc] =loglog_fit(h_values, e_implic_midpoint)
    [p_bw_loc,k_bw_loc] =loglog_fit(h_values, e_backw_e)

    figure(5)
    loglog(h_values,abs(analytical_value-X0),'MarkerFaceColor','k');
    hold on
    loglog(h_values,e_fw_e,'MarkerFaceColor','r');
    loglog(h_values,e_midpoint,'MarkerFaceColor','g');
    loglog(h_values,e_backw_e,'MarkerFaceColor','b');
    loglog(h_values,e_implic_midpoint,'MarkerFaceColor','m');

    
    hold off


end 