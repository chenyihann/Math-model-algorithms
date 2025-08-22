function eta_sb_tow2mirr_array = cal_tower_shadow_efficiency(S, every_mirr_n_unit_vector_array, x, y, z, l_array, w_array)
    R = 3.5;  
    H = 88;   
    S = S / norm(S);
    mirror_num = length(x);
    eta_sb_tow2mirr_array = ones(1, mirror_num);

    % 只计算前方镜子
    need_cal = find( [x(:), y(:), zeros(mirror_num,1)] * S' > 0 );

    % 预计算旋转矩阵
    R_all = cellfun(@(n) cal_rotation_matrix(n), num2cell(every_mirr_n_unit_vector_array,2), 'UniformOutput', false);

    x_step = 1; y_step = 1;

    for idx = need_cal'
        % 网格生成
        [Xg, Yg] = meshgrid(-l_array(idx)/2:x_step:l_array(idx)/2, -w_array(idx)/2:y_step:w_array(idx)/2);
        P_local = [Xg(:), Yg(:), zeros(numel(Xg),1)];
        P_global = (R_all{idx}' * P_local')' + [x(idx), y(idx), z(idx)];

        D = -S; % 射线方向
        a = D(1)^2 + D(2)^2;
        b = 2*(P_global(:,1)*D(1) + P_global(:,2)*D(2));
        c = P_global(:,1).^2 + P_global(:,2).^2 - R^2;
        disc = b.^2 - 4*a.*c;

        % 侧面命中
        valid_disc = disc >= 0;
        sqrt_disc = sqrt(max(disc,0));
        t1 = (-b - sqrt_disc) ./ (2*a);
        t2 = (-b + sqrt_disc) ./ (2*a);
        side_hit = false(size(disc));
        for t = [t1, t2]
            z_int = P_global(:,3) + t * D(3);
            side_hit = side_hit | (valid_disc & t >= 0 & z_int >= 0 & z_int <= H);
        end

        % 底面命中
        bottom_hit = false(size(side_hit));
        if abs(D(3)) > 1e-6
            t_bottom = -P_global(:,3) / D(3);
            x_int = P_global(:,1) + t_bottom * D(1);
            y_int = P_global(:,2) + t_bottom * D(2);
            bottom_hit = t_bottom >= 0 & (x_int.^2 + y_int.^2 <= R^2);
        end

        shaded = side_hit | bottom_hit;
        eta_sb_tow2mirr_array(idx) = 1 - nnz(shaded) / numel(shaded);
    end
end
