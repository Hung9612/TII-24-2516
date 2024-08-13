function [t_acc, t_dec, t_const, d_acc, d_dec, d_const, v_max_actual] = calculate_trajectory_parameters(v_start, v_end, v_max, a, d, x_start, x_end)
    % 

    d_total = x_end - x_start;
    d_acc = v_max^2 / (2 * a);
    d_dec = v_max^2 / (2 * d);
    d_const = d_total - d_acc - d_dec;

    % 
    if d_const < 0
        % 
        v_max_actual = sqrt((d_total * a * d) / (a + d));
        d_acc = v_max_actual^2 / (2 * a);
        d_dec = v_max_actual^2 / (2 * d);
        d_const = 0; % 
    else
        v_max_actual = v_max;
    end

    % 
    t_acc = v_max_actual / a;
    t_dec = v_max_actual / d;
    t_const = d_const / v_max_actual;

end


