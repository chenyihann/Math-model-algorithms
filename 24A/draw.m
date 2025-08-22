function draw(head_verts, body_verts, b, x_pos, y_pos)
theta = 0:0.1:32*pi;
r = b .* theta;
x = r .* cos(theta);
y = r .* sin(theta);

% 绘制螺旋线
figure;
plot(x, y ,'b-', 'LineWidth', 1.5);

hold on;

plot(x_pos,y_pos,'r.-','MarkerSize',7);  % 绘制每个把手

% 绘制龙头板凳
% patch 函数可以绘制并填充一个多边形
patch(head_verts(:,1), head_verts(:,2), 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2);

% 3. 绘制发生碰撞的龙身板凳
patch(body_verts(:,1), body_verts(:,2), 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'LineWidth', 2);


% 设置图形标题和坐标轴标签
title('螺旋线及各把手分布');
xlabel('x');
ylabel('y');
grid on;
axis equal;  % 使 x 和 y 轴比例相同，避免图像变形

hlod off;

end

