function eta_sb_onemirr = cal_eta_sb_one_mirr(need_cal_idx_array, x_array, y_array, z_array, every_mirr_n_unit_vector_array, aim_idx, w_array, l_array, S)
    % 目标镜 A 参数
    n_a = every_mirr_n_unit_vector_array(aim_idx,:);
    w_a = w_array(aim_idx);
    l_a = l_array(aim_idx);
    O_a = [x_array(aim_idx), y_array(aim_idx), z_array(aim_idx)];
    R_a = cal_rotation_matrix(n_a);

    % 一次性生成镜 A 网格点
    x_step = 1; y_step = 1;
    [Xg, Yg] = meshgrid(-l_a/2:x_step:l_a/2, -w_a/2:y_step:w_a/2);
    gridGlobal = R_a' * [Xg(:)'; Yg(:)'; zeros(1,numel(Xg))] + O_a';

    % 预计算所有旋转矩阵和中心坐标
    R_all = cellfun(@(n) cal_rotation_matrix(n), num2cell(every_mirr_n_unit_vector_array,2), 'UniformOutput', false);
    O_all = num2cell([x_array(:), y_array(:), z_array(:)], 2);

    % 遮挡检测
    N = size(gridGlobal, 2);
    shadedMask = false(1, N);
    S_col = S(:); % 列向量

    parfor k = 1:length(need_cal_idx_array)
        idx = need_cal_idx_array(k);
        Rb = R_all{idx};
        Ob = O_all{idx};
        lb = l_array(idx);
        wb = w_array(idx);

        Pb = Rb * (gridGlobal - Ob');
        Sb = Rb * S_col;

        if abs(Sb(3)) > 1e-6
            x2 = (Sb(3)*Pb(1,:) - Sb(1)*Pb(3,:)) / Sb(3);
            y2 = (Sb(3)*Pb(2,:) - Sb(2)*Pb(3,:)) / Sb(3);
            shadedMask = shadedMask | (abs(x2) <= lb/2 & abs(y2) <= wb/2);
        end
    end

    eta_sb_onemirr = 1 - nnz(shadedMask) / N;
end
