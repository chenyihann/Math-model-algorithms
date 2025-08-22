function [eta_sb_array, eta_sb_mirr_reciever_array, eta_sb_mirr_reflect_array, eta_sb_tow2mirr_array] = cal_eta_sb(month)
% 输入：
% month：输入月份
% 输出：
% eta_sb_array：得到该月份5个时间点所有定日镜的遮挡效率，5×N

% 读取整个 A 和 B 列
data = readmatrix('附件.xlsx', 'Range', 'A:B');
% 跳过第一行（标签）
data = data(2:end, :);
x = data(:,1);
y = data(:,2);
num_mirror = length(x);
z = zeros(num_mirror,1) + 4;  % 先生成全0矩阵，再全部加4


% 长宽数组6×6
w_array = zeros(num_mirror,1) + 6;
l_array = zeros(num_mirror,1) + 6;

time_array = [9, 10.5, 12, 13.5, 15];  % 单位为小时
num_time = length(time_array);

eta_sb_array = zeros(num_time, num_mirror); 
eta_sb_mirr_reciever_array = zeros(num_time,num_mirror);
eta_sb_mirr_reflect_array = zeros(num_time,num_mirror);
eta_sb_tow2mirr_array = zeros(num_time,num_mirror);


for i = 1:num_time

    time = time_array(i);
    % 调用函数，计算太阳角度
    [sin_alpha_s, cos_alpha_s, sin_beta_s,cos_beta_s, ~] = cal_alpha_and_beta(time, month);  

    % 太阳方向向量
    S_reciever = [-cos_alpha_s*sin_beta_s, -cos_alpha_s*cos_beta_s, -sin_alpha_s];

    [every_mirr_n_unit_vector_array, S_reflect_array] = cal_mirror_n(S_reciever,x,y,z);

    % 计算塔阴影对定日镜的遮挡
    eta_sb_tow2mirr_oneS_array = cal_tower_shadow_efficiency(S_reciever, every_mirr_n_unit_vector_array, x, y, z, l_array, w_array);

    eta_sb_tow2mirr_array(i,:) = eta_sb_tow2mirr_oneS_array;

    for j = 1:num_mirror  % 计算所有定日镜
        aim_idx = j;
        S_reflect = S_reflect_array(j,:);  % 形状1×3
        eta_sb_tow2mirr = eta_sb_tow2mirr_oneS_array(j);
        
        % 计算入射光的遮挡效率
        need_cal_idx_array_reciever = cal_ifcal(x,y,aim_idx,S_reciever,1,0);
        eta_sb_onemirr_reciever = cal_eta_sb_one_mirr(need_cal_idx_array_reciever,x,y,z,every_mirr_n_unit_vector_array, aim_idx,w_array,l_array,S_reciever);
        

        % 计算反射光的遮挡效率
        need_cal_idx_array_reflect = cal_ifcal(x,y,aim_idx,S_reflect,0,1);
        eta_sb_onemirr_reflect = cal_eta_sb_one_mirr(need_cal_idx_array_reflect,x,y,z,every_mirr_n_unit_vector_array, aim_idx,w_array,l_array,S_reflect);
        
        eta_sb_onemirr = eta_sb_onemirr_reciever*eta_sb_onemirr_reflect*eta_sb_tow2mirr;
        eta_sb_array(i,j) = eta_sb_onemirr;
        eta_sb_mirr_reciever_array(i,j) = eta_sb_onemirr_reciever;
        eta_sb_mirr_reflect_array(i,j) = eta_sb_onemirr_reflect;
    end
end


end


    