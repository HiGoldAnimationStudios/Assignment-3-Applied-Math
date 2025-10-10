function global_trunc_error(rate_func_in, X, t_0, t_f)

    h_values = logspace(-4, -1, 100);
    e_fw_g_e = zeros(size(h_values));
    g_e_midpoint = zeros(size(h_values));
    %disp(length(h_values))
    %disp(length(e_fw_e))

    X_0 = X(t_0);
    X_f = X(t_f);

    for i = 1:length(h_values)
        h = h_values(i);

        t_current = t_0;
        x_current_fw_e = X_0;

        while t_current < t_f
            if t_current + h > t_f
                h_step = t_f - t_current;
            else
                h_step = h;
            end

            x_current_fw_e = forward_euler_step(rate_func_in, t_current, x_current_fw_e, h_step);
            t_current = t_current + h_step;
        end

        t_current = t_0;
        x_current_mid = X_0;

        while t_current < t_f
            if t_current + h > t_f
                h_step = t_f - t_current;
            else
                h_step = h;
            end

            x_current_mid = explicit_midpoint_step(rate_func_in, t_current, x_current_mid, h_step);
            t_current = t_current + h_step;
        end

        e_fw_g_e(i) = abs(x_current_fw_e - X_f);
        g_e_midpoint(i) = abs(x_current_mid  - X_f);

    end 

    [p2_glob,k2_glob] =loglog_fit(h_values, e_fw_g_e)
    [p3_glob,k3_glob] =loglog_fit(h_values, g_e_midpoint)


    figure(67)
    hold on
    loglog(h_values,e_fw_g_e,'ro','MarkerFaceColor','r');
    loglog(h_values,g_e_midpoint,'bo','MarkerFaceColor','b');
    hold off

end 
