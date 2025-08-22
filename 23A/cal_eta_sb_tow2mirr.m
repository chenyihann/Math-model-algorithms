function eta_sb_tow2mirr_array = cal_eta_sb_tow2mirr(S,every_mirr_n_unit_vector_array,x,y,z,l_array,w_array)
% 输入：
% S：某个时间点的太阳光线方向向量size=1×3
% every_mirr_n_unit_vector_array：所有定日镜在S太阳光线下的法向量，size=N×3
% x,y,z: 定日镜中点坐标
% 输出：
% eta_sn_tow2mirr: 在特定的S下，集热器平均对定日镜的遮挡效率


R = 3.5;  % 集热器半径
H = 88;   % 总高度

S_x = S(1);
S_y = S(2);
cos_alpha = abs(S(3));  % 取锐角
alpha = acos(cos_alpha);


% 计算正交向量
if abs(S_x) < 1e-6 && abs(S_y) < 1e-6
    e1 = [1, 0, 0];
else
    e1 = [-S_y, S_x, 0];
    e1 = e1 / norm(e1);
end
e2 = [S_x, S_y, 0];
if norm(e2) < 1e-6
    e2 = [0, 1, 0];
else
    e2 = e2 / norm(e2);
end


OA_vec = R.*e1;
x1 = OA_vec(1);
y1 = OA_vec(2);

L = H*tan(alpha);

OB_vec = R.*e1 + L.*e2 + OA_vec;
x2 = OB_vec(1);
y2 = OB_vec(2);

OD_vec = -R.*e1;
x4 = OD_vec(1);
y4 = OD_vec(2);

OC_vec = 2*OD_vec + L*e2;
x3 = OC_vec(1);
y3 = OC_vec(2);


if abs(e1(1)) < 1e-6
    k = Inf; % 垂直线情况
else
    k = e1(2) / e1(1);
end

x_step = 0.5;
y_step = 0.5;

mirror_num = length(x);
d_min = 50;  % 最小距离
eta_sb_tow2mirr_array = zeros(1,mirror_num); % 预分配

for i = 1:mirror_num
    
    x0 = x(i);
    y0 = y(i);
    z0 = z(i);
    
    if isinf(k)
        d = abs(x0 - x2);
        orign_val = y2;
        aim_val = y2 - y0;
    else
        orign_val = -k * x2 + y2;
        aim_val = k * (x0 - x2) + y2 - y0;
        d = abs(k * aim_val) / sqrt(k^2 + 1);
    end

    if d <= d_min && orign_val * aim_val >= 0  % 该索引需要计算是否产生遮挡

        O_mirr_center = [x0;y0;z0];
        n0 = every_mirr_n_unit_vector_array(i,:);
        R0 = cal_rotation_matrix(n0);
        l0 = l_array(i);
        w0 = w_array(i);
        l_num = round(l0/x_step);  % 四舍五入到最接近的整数
        w_num = round(w0/y_step);
        
        N = l_num*w_num;
    
        % 使用 meshgrid 生成网格点
        x_grid = linspace(-l0/2, l0/2, l_num);
        y_grid = linspace(-w0/2, w0/2, w_num);
        [X0, Y0] = meshgrid(x_grid, y_grid);
        P0 = [X0(:), Y0(:), zeros(N, 1)]; % N×3
    
        % 转换到地面坐标系
        P1 = (R0' * P0')' + repmat(O_mirr_center', N, 1);  % N×3
    
        % 向量化阴影检查
        result = isPointInRectangle(P1(:,1), P1(:,2), x1, y1, x2, y2,x3,y3, x4, y4);
        M = sum(result);
    
        eta_sb_one_tow2mirr = 1 - M/N;

    else
        eta_sb_one_tow2mirr = 1;
    end

    eta_sb_tow2mirr_array(i) = eta_sb_one_tow2mirr;

end

end





    
   
            








