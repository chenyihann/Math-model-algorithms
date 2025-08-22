function draw_full_dragon(all_x, all_y, b)
    % 该函数在特定时刻绘制整条龙的几何形态
    
    board_width = 0.30;
    handle_dist_from_end = 0.275; % 把手(孔)距离板凳末端的距离
    
    theta = 0:0.1:32*pi;
    r = b .* theta;
    x = r .* cos(theta);
    y = r .* sin(theta);
    
    % 绘制螺旋线
    figure;
    plot(x, y ,'b-', 'LineWidth', 1.5);

    hold on;

    % 循环绘制每一节板凳 
    num_handles = length(all_x);
    for k = 1:(num_handles - 1)
        
        % 获取当前板凳的两个把手位置
        P_start = [all_x(k), all_y(k)];
        P_end   = [all_x(k+1), all_y(k+1)];
        
        % 计算板凳的单位方向向量和垂直向量 (这是正确的核心逻辑)
        vec = P_end - P_start;
        unit_vec = vec / norm(vec);
        perp_vec = [-unit_vec(2), unit_vec(1)];
        
        % 计算板凳的物理前端和后端点
        % 前端点在第一个把手P_start的前方
        bench_front_center = P_start - handle_dist_from_end * unit_vec;
        % 后端点在第二个把手P_end的后方
        bench_back_center = P_end + handle_dist_from_end * unit_vec;
        
        % 计算该板凳的四个顶点
        v1 = bench_front_center + (board_width/2) * perp_vec;
        v2 = bench_front_center - (board_width/2) * perp_vec;
        v3 = bench_back_center  - (board_width/2) * perp_vec;
        v4 = bench_back_center  + (board_width/2) * perp_vec;
        
        color = 'r'; % 红色
        
        % 使用 patch 绘制带填充的矩形
        patch([v1(1), v2(1), v3(1), v4(1)], [v1(2), v2(2), v3(2), v4(2)], ...
              color, 'FaceAlpha', 0.6, 'EdgeColor', 'k');
    end

    % 绘制把手位置 (在顶层，看得更清楚)
    plot(all_x, all_y, 'k.-', 'MarkerFaceColor', 'y', 'MarkerSize', 5);

    % --- 图形格式化 ---
    title('螺旋线和板凳龙分布');
    xlabel('x (m)');
    ylabel('y (m)');
    grid on;
    axis equal;  % 必须使用，确保x,y轴比例相同，避免几何形状失真
    hold off;
end