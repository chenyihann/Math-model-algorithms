%{
在盘入到掉头空间的边界时，龙头不能和龙身相碰。边界条件是：在掉头空间边界R=4.5处，龙头与龙身相碰。
龙头只可能与前面的板凳相撞，所以在搜索碰撞的板凳时，给到前30

螺距越小，越容易发生碰撞
%}

clc,clear;
tic;

%% 利用二分法找到最小的螺距

% 初始化二分查找边界
l = 0.4;
r = 0.55;
precision = 0.000001;

fprintf('开始二分搜索...\n');
fprintf('%-15s %-15s %-15s %-10s\n', '下界 (l)', '上界 (r)', '中点 (mid)', 'ifcrush');
fprintf('------------------------------------------------------------\n');


while (r - l) >= precision
    mid = l + (r - l) / 2;

    ifcrush = cal_crush(mid);

    fprintf('%-15.8f %-15.8f %-15.8f %-10d\n', l, r, mid, ifcrush);


    % ifcrush=2时，碰撞半径小于R，所以需要减小半径，让碰撞半径增大
    if ifcrush == 1 || ifcrush == 0  % 这个螺距过小，需要缩小搜索螺距的下限
        l = mid;
    elseif ifcrush == 2  
        r = mid;
    end
end

fprintf('------------------------------------------------------------\n');
fprintf('搜索完成。\n');
fprintf('最终找到的最小螺距为: %.5f\n', mid);

toc;


