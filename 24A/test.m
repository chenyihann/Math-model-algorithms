% 龙头板的几何参数
R = 4.5;
board_width = 0.30;
handle_dist_from_end = 0.275;
head_board_len = 3.41;

for i = 1:length(t_val)

    P0 = [x_pos(i,1), y_pos(i,1)]; % 龙头前把手
    P1 = [x_pos(i,2), y_pos(i,2)]; % 龙头后把手

    vec_P0P1 = P1-P0;
    unit_vec = vec_P0P1 / norm(vec_P0P1); % 板凳方向的单位向量
    perp_vec = [-unit_vec(2), unit_vec(1)]; % 与板凳垂直的单位向量

    % 龙头前部和后部的中心点
    center_front = P0 - handle_dist_from_end * unit_vec;
    center_back = P0 + (head_board_len - handle_dist_from_end) * unit_vec;

    % 计算龙头四个顶点
    head_v1 = center_front + (board_width/2) * perp_vec;
    head_v2 = center_front - (board_width/2) * perp_vec;
    head_v3 = center_back  - (board_width/2) * perp_vec;
    head_v4 = center_back  + (board_width/2) * perp_vec;
    head_vertices = [head_v1; head_v2; head_v3; head_v4];
    
    for j = 3:total_num

        % 计算龙身j的四个顶点
        P_j_prev = [x_pos(i, j-1), y_pos(i, j-1)];
        P_j_curr = [x_pos(i, j),   y_pos(i, j)];

        body_board_len = 2.20; % 龙身板长
        vec_body = P_j_curr - P_j_prev;
        unit_vec_body = vec_body / norm(vec_body);
        perp_vec_body = [-unit_vec_body(2), unit_vec_body(1)];

        center_front_body = P_j_prev - handle_dist_from_end * unit_vec_body;
        center_back_body  = P_j_prev + (body_board_len - handle_dist_from_end) * unit_vec_body;
        
        body_v1 = center_front_body + (board_width/2) * perp_vec_body;
        body_v2 = center_front_body - (board_width/2) * perp_vec_body;
        body_v3 = center_back_body  - (board_width/2) * perp_vec_body;
        body_v4 = center_back_body  + (board_width/2) * perp_vec_body;
        body_vertices = [body_v1; body_v2; body_v3; body_v4];

        if check_collision(head_vertices, body_vertices)
            disp(['碰撞的时间为：', num2str(i*time_step+410)]);
            theta0 = theta0_values(i);
            r0 = b * theta0;
            draw(b,x_pos(i,:),y_pos(i,:));  % 调用绘图函数, 当前时间
            if r0 >= R
                disp(['碰撞位置大于R，此时的螺距为：', num2str(d)]);
                disp(['碰撞半径r为：', num2str(r0)]);
                ifcrush = 1;
                return;
            else
                disp(['碰撞位置小于R，此时的螺距为：', num2str(d)])
                disp(['碰撞半径r为：', num2str(r0)]);
                ifcrush = 2;
                return;
            end
        end
    end
end


%% 判断在R=4.5时是否发生碰撞
R = 4.5;
d_lim = 0.15;
Q1P0 = sqrt((0.15)^2 + (0.275)^2);
P1P0 = 2.86;
a = Q1P0/P1P0;
M_cos = a*(0.275 / Q1P0);
M_sin = a*(0.15/Q1P0);
M = [M_cos, -M_sin; M_sin, M_cos];

x_q1_pos = zeros(1,length(t_val));
y_q1_pos = zeros(1,length(t_val));
x_q2_pos = zeros(1,length(t_val));
y_q2_pos = zeros(1,length(t_val));

for i = 1:length(t_val)
    % 第0个把手的坐标
    x0 = x_pos(i,1);
    y0 = y_pos(i,1);
    % 第1个把手的坐标
    x1 = x_pos(i,2);
    y1 = y_pos(i,2);
    delta = [x0-x1;y0-y1];

    % 计算龙头一个顶点的坐标Q1

    result1 = M'*delta + [x0;y0];  % M要转置
    x_q1 = result1(1);
    y_q1 = result1(2);
    x_q1_pos(i) = x_q1;  % 储存进数组
    y_q1_pos(i) = y_q1;

    % 计算Q2坐标
    result2 = M*delta + [x0;y0];
    x_q2 = result2(1);
    y_q2 = result2(2);
    x_q2_pos(i) = x_q2;
    y_q2_pos(i) = y_q2;
end


for i = 1:length(t_val)
    % Q1和Q2的坐标
    x_q1 = x_q1_pos(i);
    y_q1 = y_q1_pos(i);
    x_q2 = x_q2_pos(i);
    y_q2 = y_q2_pos(i);
    for j = 3:total_num
        % 龙头不会与自身相撞，所以跳过j=1和2的情况，从第三个把手开始计算
        k = k_list(i,j);
        x_k = x_pos(i,j);
        y_k = y_pos(i,j);
        d1 = abs(k*(x_q1 - x_k) + y_k - y_q1) / sqrt(k^2 + 1);
        d2 = abs(k*(x_q2 - x_k) + y_k - y_q2) / sqrt(k^2 + 1);
        if d1 <= d_lim || d2 <= d_lim
            disp(['碰撞的时间为：', num2str(i*time_step)]);
            theta0 = theta0_values(i);  % 龙头第一个把手的角度theta0
            r0 = b * theta0;  % 第一个把手的半径
            draw(b, x_pos(i,:),y_pos(i,:));
            if r0 >= R
                ifcrush = 1;
                disp(['此时的螺距为：', num2str(b*2*pi)]);
                disp(['此时的半径为：', num2str(r0)]);
                return;
            elseif r0 < R
                ifcrush = 2;
                disp(['此时的螺距为：', num2str(b*2*pi)]);
                disp(['此时的半径为：', num2str(r0)]);
                return;
            end
        end
    end
end
