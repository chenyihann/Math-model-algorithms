function [table1, table2] = q1zhu()
    % 定义常量
    phi = 39.4 * pi / 180;  % 纬度（弧度）
    H_alt = 3;              % 海拔（km）
    G0 = 1.366;             % 太阳常数（kW/m²）
    h_tower = 80;           % 吸收塔高度（m）
    h_mirror = 4;           % 定日镜安装高度（m）
    mirror_area = 36;        % 单面定日镜面积（m²）
    center = [0;0;84];
    
    % 计算DNI参数
    a = 0.4237 - 0.00821 * (6 - H_alt)^2;
    b = 0.5055 + 0.00595 * (6.5 - H_alt)^2;
    c = 0.2711 + 0.01858 * (2.5 - H_alt)^2;
    
    % 读取定日镜位置
    data = xlsread('附件.xlsx', 'Sheet1');
    x = data(:, 1);  % x坐标(m)
    y = data(:, 2);  % y坐标(m)
    z = zeros(length(x),1) + 4;
    X = [x,y,z];
    w_array = zeros(length(x),1) + 6;  % 宽度数组
    l_array = zeros(length(x),1) + 6;  % 长度数组
    N = length(x);                     % 定日镜数量
    total_area = N * mirror_area;      % 总镜面面积
    
    % 计算每面定日镜到集热器的距离
    d_HR = sqrt(x.^2 + y.^2 + (h_tower - h_mirror)^2);
    % 计算大气透射率
    eta_at = 0.99321 - 0.0001176 * d_HR + 1.97e-8 * d_HR.^2;
    
    % 月份设置（每月21日的年积日）
    doy_list = [21, 52, 80, 111, 141, 172, 202, 233, 264, 294, 325, 355];
    D = mod(doy_list - 80, 365);  % 以春分（3月21日）为基准
    times = [9, 10.5, 12, 13.5, 15]; % 当地时间（小时）
    
    % 初始化结果存储
    monthly_avg_eta_opt = zeros(12, 1);
    monthly_avg_eta_cos = zeros(12, 1);
    monthly_avg_eta_sb = ones(12, 1);  % 阴影遮挡效率假设为1
    monthly_avg_eta_trunc = ones(12, 1); % 截断效率假设为1
    monthly_avg_E_field = zeros(12, 1);
    monthly_avg_E_per_area = zeros(12, 1);
    monthly_avg_eta_sb_mirr_reciever_array = zeros(12,1);
    monthly_avg_eta_sb_mirr_reflect_array = zeros(12,1);
    monthly_avg_eta_sb_tow2mirr_array = zeros(12,1);
    
    % 主循环：遍历每个月
    for i_month = 1:12
        D_i = D(i_month);
        % 初始化该月各时刻的数据
        eta_opt_time = zeros(5, 1);
        eta_cos_time = zeros(5, 1);
        E_field_time = zeros(5, 1);
        eta_truc_time = zeros(5, 1);


        [eta_sb_array, eta_sb_mirr_reciever_array, eta_sb_mirr_reflect_array, eta_sb_tow2mirr_array] = cal_eta_sb(i_month);  

        
        % 遍历每天5个时刻
        for i_time = 1:5
            ST = times(i_time);
            omega = (pi / 12) * (ST - 12); % 太阳时角
            
            % 计算太阳赤纬角
            delta = asin(sin(2 * pi * D_i / 365) * sin(deg2rad(23.45)));
            % 计算太阳高度角
            sin_alpha_s = cos(delta) * cos(phi) * cos(omega) + sin(delta) * sin(phi);
            
            if sin_alpha_s <= 0
                % 太阳在地平线下，无输出
                eta_opt_time(i_time) = 0;
                eta_cos_time(i_time) = 0;
                E_field_time(i_time) = 0;
                continue;
            end
            
            alpha_s = asin(sin_alpha_s);
            % 计算太阳方位角
            cos_gamma_s = (sin(delta) - sin(alpha_s) * sin(phi)) / (cos(alpha_s) * cos(phi));
            cos_gamma_s = max(min(cos_gamma_s, 1), -1); % 确保在[-1,1]范围内
            gamma_s = acos(cos_gamma_s);
            if omega >= 0
                gamma_s = 2 * pi - gamma_s; % 下午调整方位角
            end
            
            % 计算DNI
            DNI = G0 * (a + b * exp(-c / sin(alpha_s)));
            if DNI < 0
                DNI = 0;
            end
            
            % 计算太阳方向向量
            S = [-cos(alpha_s) * sin(gamma_s), -cos(alpha_s) * cos(gamma_s), -sin(alpha_s)];
            % 初始化该时刻定日镜数据
            eta_opt_mirrors = zeros(N, 1);
            eta_cos_mirrors = zeros(N, 1);
            % 该时刻下的截断效率
            Y = jieduan(X,l_array,w_array,alpha_s,gamma_s,center);
            
            % 遍历每面定日镜
            for i_mirror = 1:N
                % 反射方向向量（归一化）
                R_vec = [-x(i_mirror), -y(i_mirror), h_tower - h_mirror];
                R = R_vec / norm(R_vec);
                
                % 计算法向量
                sum_vec = R - S;
                norm_sum = norm(sum_vec);
                if norm_sum < 1e-6
                    eta_cos = 0; % 法向量无效
                else
                    N_vec = sum_vec / norm_sum;
                    cos_theta = dot(S, N_vec);
                    eta_cos = max(cos_theta, 0); % 确保非负
                end
                
                % 计算光学效率
                eta_trunc = Y(i_mirror);
                eta_ref = 0.92;
                eta_sb = eta_sb_array(i_time, i_mirror);
                eta_opt = eta_cos * eta_at(i_mirror) * eta_ref * eta_trunc * eta_sb; % η_ref = 0.92
                eta_opt_mirrors(i_mirror) = eta_opt;
                eta_cos_mirrors(i_mirror) = eta_cos;
            end
            
            % 该时刻镜场总输出热功率（kW）
            E_field_time(i_time) = DNI * mirror_area * sum(eta_opt_mirrors);
            % 该时刻平均光学效率
            eta_opt_time(i_time) = mean(eta_opt_mirrors);
            eta_cos_time(i_time) = mean(eta_cos_mirrors);
            eta_truc_time(i_time) = mean(Y);
        end
        
        % 储存每个月的数据
        monthly_avg_eta_trunc = mean(eta_truc_time);
        monthly_avg_eta_sb(i_month) = mean(eta_sb_array(:));
        monthly_avg_eta_sb_mirr_reciever_array(i_month) = mean(eta_sb_mirr_reciever_array(:));
        monthly_avg_eta_sb_mirr_reflect_array(i_month) = mean(eta_sb_mirr_reflect_array(:));
        monthly_avg_eta_sb_tow2mirr_array(i_month) = mean(eta_sb_tow2mirr_array(:));
        monthly_avg_eta_opt(i_month) = mean(eta_opt_time);
        monthly_avg_eta_cos(i_month) = mean(eta_cos_time);
        monthly_avg_E_field(i_month) = mean(E_field_time); % kW
        monthly_avg_E_per_area(i_month) = monthly_avg_E_field(i_month) / total_area;
    
    end
    
    % 计算年平均指标
    annual_avg_eta_sb = mean(monthly_avg_eta_sb(:));
    annual_avg_eta_opt = mean(monthly_avg_eta_opt);
    annual_avg_eta_cos = mean(monthly_avg_eta_cos);
    annual_avg_eta_trunc = mean(monthly_avg_eta_trunc);
    annual_avg_E_field = mean(monthly_avg_E_field) / 1000; % 转换为MW
    annual_avg_E_per_area = annual_avg_E_field * 1000 / total_area; % kW/m²


    
    % 构建表1（月平均数据）
    months = {'1月21日'; '2月21日'; '3月21日'; '4月21日'; '5月21日'; '6月21日'; ...
              '7月21日'; '8月21日'; '9月21日'; '10月21日'; '11月21日'; '12月21日'};
    table1 = table(months, monthly_avg_eta_opt, monthly_avg_eta_cos, ...
                   monthly_avg_eta_sb, monthly_avg_eta_trunc, ...
                   monthly_avg_E_per_area, ...
                  'VariableNames', {'Date', 'AvgOpticalEfficiency', ...
                  'AvgCosineEfficiency', 'AvgShadingEfficiency', ...
                  'AvgTruncationEfficiency', 'AvgOutputPerArea'});
    
    % 构建表2（年平均数据）
    annual_data = {annual_avg_eta_opt, annual_avg_eta_cos, annual_avg_eta_sb, ...
                   annual_avg_eta_trunc, annual_avg_E_field, annual_avg_E_per_area};
    table2 = table(annual_data{1}, annual_data{2}, annual_data{3}, annual_data{4}, ...
                  annual_data{5}, annual_data{6}, ...
                  'VariableNames', {'AnnualAvgOpticalEfficiency', ...
                  'AnnualAvgCosineEfficiency', 'AnnualAvgShadingEfficiency', ...
                  'AnnualAvgTruncationEfficiency', 'AnnualAvgOutputPower_MW', ...
                  'AnnualAvgOutputPerArea_kW_per_m2'});
    
    % 显示结果
    disp('表1 问题1每月21日平均光学效率及输出功率');
    disp(table1);
    disp('表2 问题1年平均光学效率及输出功率');
    disp(table2);
end