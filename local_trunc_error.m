function local_trunc_error(rate_func_in, X, t_ref)
    h_values = logspace(-5,1,100);
    
    e_fw_e = zeros(size(h_values));
    
    e_midpoint = zeros(size(h_values));
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

    end
    
    disp(size(h_values))
    disp(size(analytical_value))
    disp(size(e_fw_e))
    disp(size(e_midpoint))

 
    [p1_loc,k1_loc] = loglog_fit(h_values, abs(analytical_value-X0))
    [p2_loc,k2_loc] =loglog_fit(h_values, e_fw_e)
    [p3_loc,k3_loc] =loglog_fit(h_values, e_midpoint)

    figure(5)
    loglog(h_values,abs(analytical_value-X0),'go','MarkerFaceColor','g');
    hold on
    loglog(h_values,e_fw_e,'ro','MarkerFaceColor','r');
    loglog(h_values,e_midpoint,'bo','MarkerFaceColor','b');

    
    hold off


end 