%{ 
判断龙头的两个顶点与板中线的距离，如果小于d0则判定为碰撞。如何说明龙头后两个顶点不会与龙身发生碰撞？
按照Q1的思路，不断的增大t，当d1和d2有小于d0时，记录时间t。初始步长取为1，找到大致时间区间再逐步缩小。
步长缩小后，还需要从头迭代，置时间满足精度要求。
可以只计算前20个把手，只可能与前面的板凳发生碰撞，减少运行时间。（说明）

同心圆近似模判断碰撞的时间区间

以下代码不采用同心圆模拟，在一定时间范围内搜索。可以使用二分法加快搜索。
%} 

clc,clear;

total_num = 20;  % 把手数目

options = optimoptions('fsolve', 'Display','off', 'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9);

tic;
%% 求解theta与t的关系式

% 初始化全局变量
b = 0.55/(2*pi);
time_lower = 400;
time_upper = 415;
t_val0 = 0:0.1:time_lower;
time_step = 0.001;  % 搜索时间步长
t_val = time_lower:time_step:time_upper;    % 搜索时间区间
C = -b*0.5*(32*pi*sqrt((32*pi)^2 + 1) + log(32*pi + sqrt((32*pi)^2 + 1)));  % 积分常数

% 使用~忽略输出
[theta0_initial, ~] = cal_theta0(t_val0, 32*pi, b, C, options);  % 计算得到第410s处的theta0
[~, theta0_values] = cal_theta0(t_val, theta0_initial, b, C, options);  % 将410-413时间区间带入计算得到对应的theta0

theta_sol = zeros(length(t_val), total_num); % 初始化数组，存储所有解

%% 求解各个把手theta的递推关系式

for i = 1:length(t_val)
    for j = 1:total_num
        if j == 1
            theta_val = theta0_values(i);
        else
            if j == 2
                L = 2.86;
            else
                L = 1.65;
            end
            
            theta_prev = theta_sol(i, j-1);  % 使用同一个时间，前一个把手的解作为初值

            r_prev = b*theta_prev;

            theta_guess = theta_prev + L/r_prev;
            
            fun2 = @(theta_k) b^2*(theta_prev^2 + theta_k^2) - L^2 - 2*b^2*theta_k*theta_prev*cos(theta_k - theta_prev);
            theta_val = fsolve(fun2, theta_guess, options);
        end 
        theta_sol(i, j) = theta_val; % 索引从1开始
    end
end

%% 求解所有把手在410-413s时的位置
x_pos = b*theta_sol.*cos(theta_sol);
y_pos = b*theta_sol.*sin(theta_sol);

%% 速度的求解
v_sol = zeros(length(t_val), total_num); % 初始化数组，存储所有解
k_list = zeros(length(t_val), total_num);
for i = 1:length(t_val)
    for j = 1:total_num
        if j==1
            v_val = 1;
            v_sol(i,j) = v_val;
        else
            theta1 = theta_sol(i,j-1);
            theta2 = theta_sol(i,j);
            k1 = (b*sin(theta1) + b*theta1*cos(theta1))/(b*cos(theta1)-b*theta1*sin(theta1));
            k2 = (b*sin(theta2) + b*theta2*cos(theta2))/(b*cos(theta2)-b*theta2*sin(theta2));
            k = (y_pos(i,j)-y_pos(i,j-1))/(x_pos(i,j)-x_pos(i,j-1));
            k_list(i,j) = k;
            apha = atan(abs((k1-k) / (1 + k1*k)));
            beta = atan(abs((k2-k) / (1 + k2*k)));  
            v_prev = v_sol(i,j-1);
            fun3 = @(v) v*cos(beta) - v_prev*cos(apha);
            v_val = fsolve(fun3,v_prev, options); 
            v_sol(i,j) = v_val;
        end
    end
end
         
%% 判断是否发生碰撞

board_width = 0.30;
handle_dist_from_end = 0.275;
% head_board_len = 3.41;

for i = 1:length(t_val)

    % 创建标志变量
    flag = 0;

    P0 = [x_pos(i,1), y_pos(i,1)]; % 龙头前把手
    P1 = [x_pos(i,2), y_pos(i,2)]; % 龙头后把手

    vec_P0P1 = P1-P0;
    unit_vec = vec_P0P1 / norm(vec_P0P1); % 板凳方向的单位向量
    perp_vec = [-unit_vec(2), unit_vec(1)]; % 与板凳垂直的单位向量

    % 龙头前部和后部的中心点
    center_front = P0 - handle_dist_from_end * unit_vec;
    center_back = P1 + handle_dist_from_end * unit_vec;

    % 计算龙头四个顶点
    head_v1 = center_front + (board_width/2) * perp_vec;
    head_v2 = center_front - (board_width/2) * perp_vec;
    head_v3 = center_back  - (board_width/2) * perp_vec;
    head_v4 = center_back  + (board_width/2) * perp_vec;
    head_vertices = [head_v1; head_v2; head_v3; head_v4];
    
    for j = 4:total_num  % 从第2节龙身开始计算，因为龙头一定与第1节龙身有重叠

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
            disp(['碰撞的时间为：', num2str(i*time_step+time_lower)]);
            theta0 = theta0_values(i);
            r0 = b * theta0;
            disp(['碰撞半径：', num2str(r0)]);
            draw_full_dragon(x_pos(i,:), y_pos(i,:), b);  % 调用绘图函数, 当前时间
            flag = 1;
            break;
        end
    end
    if flag == 1
        break;
    end
end

toc;
