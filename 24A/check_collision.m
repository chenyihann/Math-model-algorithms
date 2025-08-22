function is_colliding = check_collision(rect1_vertices, rect2_vertices)
    % rect1_vertices 和 rect2_vertices 都是 4x2 的矩阵，每行是一个顶点的[x, y]坐标
    
    % 获取两个矩形的四条边的法向量作为分离轴
    axes = zeros(4, 2);
    
    % 矩形1的轴
    edge1 = rect1_vertices(2,:) - rect1_vertices(1,:);
    axes(1,:) = [-edge1(2), edge1(1)]; % 法向量
    edge2 = rect1_vertices(3,:) - rect1_vertices(2,:);
    axes(2,:) = [-edge2(2), edge2(1)]; % 法向量
    
    % 矩形2的轴
    edge3 = rect2_vertices(2,:) - rect2_vertices(1,:);
    axes(3,:) = [-edge3(2), edge3(1)]; % 法向量
    edge4 = rect2_vertices(3,:) - rect2_vertices(2,:);
    axes(4,:) = [-edge4(2), edge4(1)]; % 法向量

    for i = 1:4
        axis = axes(i,:);
        axis = axis / norm(axis); % 单位化
        
        % 将两个矩形的顶点投影到轴上
        proj1 = rect1_vertices * axis';
        proj2 = rect2_vertices * axis';
        
        min1 = min(proj1); max1 = max(proj1);
        min2 = min(proj2); max2 = max(proj2);
        
        % 检查投影是否分离
        if max1 < min2 || max2 < min1
            % 找到了分离轴，说明没有碰撞
            is_colliding = false;
            return;
        end
    end
    
    % 如果所有轴都无法分离它们，说明发生了碰撞
    is_colliding = true;
end