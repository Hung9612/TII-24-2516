function trajectory_points = generate_trajectory_points(v_start, v_end, v_max, a, d, x_start, x_end)
% 目标点一定要比初始点的数值要大。
% 计算轨迹参数
    [t_acc, t_dec, t_const, d_acc, d_dec, d_const, v_max_actual] = calculate_trajectory_parameters(v_start, v_end, v_max, a, d, x_start, x_end);

    % 初始化轨迹点数组
    trajectory_points = [];
    delta_t=0.001;
    % 加速阶段
    for t = 0: delta_t: t_acc-delta_t
        v = a * t;
        x = x_start + 0.5 * a * t^2;
        trajectory_points = [trajectory_points; x, v];
    end

    % 匀速阶段
    if t_const > 0
        x_lastc = x_start + d_acc;
        for t = t_acc: delta_t: t_const+t_acc-delta_t
            v = v_max_actual;
            x = x_lastc + v * (t - t_acc);
            trajectory_points = [trajectory_points; x, v];
        end
    end

    % 减速阶段
    x_last = x_start + d_acc + d_const;
    for t = t_const+t_acc: delta_t: t_dec+t_const+t_acc - delta_t
        v = v_max_actual - d * (t-(t_const+t_acc));
        x = x_last + v_max_actual*(t-(t_const+t_acc)) - 0.5 * d * (t-(t_const+t_acc))^2;
        trajectory_points = [trajectory_points; x, v];
    end

end