function [q_P2P,dq_P2P] = trajectory_points(x_start, x_end, v_start, v_end, v_max, acc, des)

    % 计算轨迹参数
    [t_acc, t_dec, t_const, d_acc, d_dec, d_const, v_max_actual] = calculate_trajectory_parameters(v_start, v_end, v_max, acc, des, x_start, x_end);

    % 初始化轨迹点数组
    Num_acc = round(t_acc*1000);
    Num_const = round(t_const*1000);
    Num_des = round(t_dec*1000);
    N= Num_acc+Num_const+Num_des;
    delta_t=0.001;
    trajectory_points = zeros(7000,2);
    % 加速阶段
    i=1;
    for t = 0: delta_t: t_acc-delta_t
        v = acc * t;
        x = x_start + 0.5 * acc * t^2;
        trajectory_points(i,:) = [x, v];
        i=i+1;
    end

    % 匀速阶段
    if t_const > 0
        x_lastc = x_start + d_acc;
        for t = t_acc: delta_t: t_const+t_acc-delta_t
            v = v_max_actual;
            x = x_lastc + v * (t - t_acc); %在加速阶段的位移上增加
            trajectory_points(i,:) = [x, v];
            i=i+1;
        end
    end

    % 减速阶段
    x_last = x_start + d_acc + d_const;
    for t = t_const+t_acc: delta_t: t_dec+t_const+t_acc - delta_t
        v = v_max_actual - des * (t-(t_const+t_acc));
        
        x = x_last + v_max_actual*(t-(t_const+t_acc))...
            - 0.5 * des * (t-(t_const+t_acc))^2; %在加速-匀速两个阶段的位移上增加
        trajectory_points(i,:) = [x, v];
        i=i+1;
    end
    for k = i:1: 7000 - N
        trajectory_points(k,:) = trajectory_points(i-1,:);
    end
q_P2P = trajectory_points(:,1);
dq_P2P = trajectory_points(:,2);
end