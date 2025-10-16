function local_trunc_error(rate_func_in, X, t_ref)
    h_values = logspace(-5,1,100);
    
    e_fw_e = zeros(size(h_values));
    e_midpoint = zeros(size(h_values));
    e_implicit_mid = zeros(size(h_values));
    e_backw_mid = zeros(size(h_values));

    analytical_value = zeros(size(h_values));

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
    
    [p_x_loc,k_x_loc] = loglog_fit(h_values, abs(analytical_value-X0))
    [p_fw_loc,k_fw_loc] =loglog_fit(h_values, e_fw_e)
    [p_ex_loc,k_ex_loc] =loglog_fit(h_values, e_midpoint)
    %[p_im_loc,k_im_loc] =loglog_fit(h_values, e_implic_midpoint)
    %[p_bw_loc,k_bw_loc] =loglog_fit(h_values, e_backw_e)

    hh = logspace(log10(min(h_values)), log10(max(h_values)), 200);



    figure() 
    loglog(h_values,abs(analytical_value-X0), 'o', 'MarkerSize',5, 'Color',[0 0 0], 'DisplayName','|X(t{+}h)-X(t)| (exact)'); 
    hold on 
    loglog(h_values,e_fw_e,'o', 'MarkerSize',5, 'Color',[1 0 0], 'DisplayName','Forward Euler Error'); 
    loglog(h_values,e_midpoint,'o', 'MarkerSize',5, 'Color',[0 1 0], 'DisplayName','Explicit Midpoint Error'); 
    
    loglog(hh, k_x_loc.*hh.^p_x_loc, 'Color','k', 'LineWidth',1.5, 'DisplayName', sprintf('Exact fit: p=%.2f', p_x_loc));
    loglog(hh, k_fw_loc .*hh.^p_fw_loc, 'Color','r', 'LineWidth',1.5, 'DisplayName', sprintf('Forward Euler fit: p=%.2f', p_fw_loc));
    loglog(hh, k_ex_loc .*hh.^p_ex_loc, 'Color','g', 'LineWidth',1.5, 'DisplayName', sprintf('Explicit Midpoint fit: p=%.2f', p_ex_loc));
    
    xlabel('step size h');
    ylabel('local error at t_{ref}{+}h');
    title(sprintf('Local Truncation Error (explicit & exact) vs. h (t_{ref} = %.3g)', t_ref));
    legend('Location','best');
    hold off

    figure() 
    loglog(h_values,e_fw_e, 'o','Color','r','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Forward Euler Error'); 
    hold on
    loglog(h_values,e_midpoint, 'o', 'Color','g','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Explicit Midpoint Error'); 
    loglog(h_values,e_backw_e, 'o', 'Color','b','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Backward Euler Error'); 
    loglog(h_values,e_implic_midpoint, 'o', 'Color','y','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Implicit Midpoint Error'); 
    xlabel('step size h');
    ylabel('local error at t_{ref}{+}h');
    title(sprintf('Local Truncation Error (implicit & explicit) vs. h (t_{ref} = %.3g)', t_ref));
    legend('Location','best');
    hold off
end 