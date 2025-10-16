function global_trunc_error(rate_func_in, X, t_0, t_f)

    h_values = logspace(-4, -1, 50);
    g_e_e_fw = zeros(size(h_values));
    g_e_ex_m = zeros(size(h_values));
    g_e_e_bw = zeros(size(h_values));
    g_e_im_m = zeros(size(h_values));
    

    X_0 = X(t_0);
    X_f = X(t_f);
    num_evals_fw_list=zeros(size(h_values));
    num_evals_ex_list=zeros(size(h_values));
    num_evals_bw_list=zeros(size(h_values));
    num_evals_im_list=zeros(size(h_values));

    for i = 1:length(h_values)
        h = h_values(i);

        t_current = t_0;
        x_current_fw_e = X_0;
        x_current_ex_mid = X_0;
        x_current_bw_e = X_0;
        x_current_im_mid = X_0;

        num_evals_fw_sum=0;
        num_evals_ex_sum=0;
        num_evals_bw_sum=0;
        num_evals_im_sum=0;

        while t_current < t_f
            if t_current + h > t_f
                h_step = t_f - t_current;
            else
                h_step = h;
            end

            [x_current_fw_e,num_evals_fw] = forward_euler_step(rate_func_in, t_current, x_current_fw_e, h_step);
            [x_current_ex_mid,num_evals_ex] = explicit_midpoint_step(rate_func_in, t_current, x_current_ex_mid, h_step);
            [x_current_bw_e,num_evals_bw] = backward_euler_step(rate_func_in,t_current,x_current_bw_e, h_step);
            [x_current_im_mid,num_evals_im] = implicit_midpoint_step(rate_func_in, t_current, x_current_im_mid, h_step);


            num_evals_fw_sum=num_evals_fw_sum+num_evals_fw;
            num_evals_ex_sum=num_evals_ex_sum+num_evals_ex;
            num_evals_bw_sum=num_evals_bw_sum+num_evals_bw;
            num_evals_im_sum=num_evals_im_sum+num_evals_im;
            t_current = t_current + h_step;
        end

        g_e_e_fw(i) = abs(x_current_fw_e - X_f);
        g_e_ex_m(i) = abs(x_current_ex_mid  - X_f);
        g_e_e_bw(i) = abs(x_current_bw_e - X_f);
        g_e_im_m(i) = abs(x_current_im_mid  - X_f);

        num_evals_fw_list(i)=num_evals_fw_sum;
        num_evals_ex_list(i)=num_evals_ex_sum;
        num_evals_bw_list(i)=num_evals_bw_sum;
        num_evals_im_list(i)=num_evals_im_sum;
    end 


    hh   = logspace(log10(min(h_values)), log10(max(h_values)), 200);

    [pfw_glob,kfw_glob] =loglog_fit(h_values, g_e_e_fw)
    [pex_glob,kex_glob] =loglog_fit(h_values, g_e_ex_m)
    [pbw_glob,kbw_glob] =loglog_fit(h_values, g_e_e_bw)
    [pim_glob,kim_glob] =loglog_fit(h_values, g_e_im_m)

    figure()
    loglog(h_values,g_e_e_fw,'o', 'MarkerSize',5, 'Color',[1 0 0], 'DisplayName','Forward Euler Error');
    hold on
    loglog(h_values,g_e_ex_m,'o', 'MarkerSize',5, 'Color',[0 1 0], 'DisplayName','Explicit Midpoint Error');

    loglog(hh, kfw_glob .*hh.^pfw_glob, 'Color','r', 'LineWidth',1.5, 'DisplayName', sprintf('Forward Euler fit: p=%.2f', pfw_glob));
    loglog(hh, kex_glob .*hh.^pex_glob, 'Color','g', 'LineWidth',1.5, 'DisplayName', sprintf('Explicit Midpoint fit: p=%.2f', pex_glob));
    xlabel('step size h');
    ylabel('global error at t_{0}{+}h');
    title(sprintf('Global Truncation Error (explicit with fit lines) vs. h (t_{f} = %.3g)', t_f));
    legend('Location','best');
    hold off

    figure()
    loglog(h_values,g_e_e_fw,'o','Color','r','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Forward Euler Error');
    hold on
    loglog(h_values,g_e_ex_m,'o', 'Color','g','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Explicit Midpoint Error');
    loglog(h_values,g_e_e_bw,'o', 'Color','b','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Backward Euler Error');
    loglog(h_values,g_e_im_m,'o', 'Color','y','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Implicit Midpoint Error');    
    xlabel('step size h');
    ylabel('global error at t_{0}{+}h');
    title(sprintf('Global Truncation Error (explicit & implicit) vs. h (t_{f} = %.3g)', t_f));
    legend('Location','best');
    hold off




    disp("with num_evals")
    
    [pfw_glob,kfw_glob] =loglog_fit(num_evals_fw_list, g_e_e_fw)
    [pex_glob,kex_glob] =loglog_fit(num_evals_ex_list, g_e_ex_m)
    [pbw_glob,kbw_glob] =loglog_fit(num_evals_bw_list, g_e_e_bw)
    [pim_glob,kim_glob] =loglog_fit(num_evals_im_list, g_e_im_m)

    
    
    
    figure()
    loglog(num_evals_fw_list,g_e_e_fw,'o', 'MarkerSize',5, 'Color',[1 0 0], 'DisplayName','Forward Euler Error');
    hold on
    loglog(num_evals_ex_list,g_e_ex_m,'o', 'MarkerSize',5, 'Color',[0 1 0], 'DisplayName','Explicit Midpoint Error');

    loglog(hh, kfw_glob .*hh.^pfw_glob, 'Color','r', 'LineWidth',1.5, 'DisplayName', sprintf('Forward Euler fit: p=%.2f', pfw_glob));
    loglog(hh, kex_glob .*hh.^pex_glob, 'Color','g', 'LineWidth',1.5, 'DisplayName', sprintf('Explicit Midpoint fit: p=%.2f', pex_glob));
    xlabel('step size h');
    ylabel('global error at t_{0}{+}h');
    title(sprintf('Global Truncation Error (explicit with fit lines) vs. number of function calls (t_{f} = %.3g)', t_f));
    legend('Location','best');
    hold off

    figure()
    loglog(num_evals_fw_list,g_e_e_fw,'o','Color','r','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Forward Euler Error');
    hold on
    loglog(num_evals_ex_list,g_e_ex_m,'o', 'Color','g','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Explicit Midpoint Error');
    loglog(num_evals_bw_list,g_e_e_bw,'o', 'Color','b','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Backward Euler Error');
    loglog(num_evals_im_list,g_e_im_m,'o', 'Color','y','MarkerSize',5, 'LineWidth',0.5, 'DisplayName','Implicit Midpoint Error');    
    xlabel('number of function calls');
    ylabel('global error at t_{0}{+}h');
    title(sprintf('Global Truncation Error (explicit & implicit) vs. number of error calls (t_{f} = %.3g)', t_f));
    legend('Location','best');
    hold off



end 
