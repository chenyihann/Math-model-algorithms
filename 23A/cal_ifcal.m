% 子函数：判断哪两个定日镜需要两两计算，距离小于
function need_cal_idx_array= cal_ifcal(x,y,i,S,reciever, reflect)
    % 输出：0或1，是否需要两两计算
    % 输入：所有定日镜的坐标，当前计算的定日镜的索引i
    % 判断是否和太阳方向同向，如果向量AB的xy坐标和太阳光线同向，且AB的模长小于50，则需要计算是否遮挡

    %当前定日镜坐标
    x0 = x(i);
    y0 = y(i);
    % 初始化
    need_cal_idx_array = [];
    
    for j = 1:length(x)
        if i ~= j
            x1 = x(j);
            y1 = y(j);
            d = sqrt((x0-x1)^2 + (y0-y1)^2);
            if reciever == 1
                vec = [x0-x1,y0-y1];
            end
            if reflect == 1
                vec = [x1-x0,y1-y0];
            end
            if (vec(1)*S(1)>0 && vec(2)*S(2)>0 ) && d <= 50
                need_cal_idx_array = [need_cal_idx_array, j];  % 将索引j加入进数组
            end
        end
    end
end
